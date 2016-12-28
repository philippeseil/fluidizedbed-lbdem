/*
 * This file is part of the LBDEMcoupling software.
 *
 * LBDEMcoupling is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright 2014 Johannes Kepler University Linz
 *
 * Author: Philippe Seil (philippe.seil@jku.at)
 */

#define MRT_USE_TRT_RELAXATION

#include "palabos3D.h"
#include "palabos3D.hh"

#include "plb_ib.h"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>

#include "inputParser.h"

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::MRTD3Q19Descriptor
#define BASEDYNAMICS MRTdynamics<T, DESCRIPTOR>( parameters.getOmega() )
#define SMAGO_BASEDYNAMICS SmagorinskyMRTdynamics<T, DESCRIPTOR>( parameters.getOmega(), cSmago )

#define DYNAMICS IBcompositeDynamics<T, DESCRIPTOR>( new BASEDYNAMICS, parameters.getOmega() )
#define SMAGO_DYNAMICS IBcompositeDynamics<T,DESCRIPTOR>( new SMAGO_BASEDYNAMICS, parameters.getOmega() )

//#define DYNAMICS BGKdynamics<T,DESCRIPTOR>(parameters.getOmega())
void writeVTK(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
              IncomprFlowParam<T> const& parameters,
              PhysUnits3D<T> const& units, plint iter)
{
  
  T p_fact = units.getPhysForce(1)/pow(units.getPhysLength(1),2)/3.;
  
  std::string fname(createFileName("vtk", iter, 6));
  
  VtkImageOutput3D<T> vtkOut(fname, units.getPhysLength(1));
  vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", units.getPhysVel(1));  
  
  MultiScalarField3D<T> p(*computeDensity(lattice));
  subtractInPlace(p,1.);
  vtkOut.writeData<float>(p,"pressure",p_fact ); 
  
  pcout << "wrote " << fname << std::endl;
}

void writePressure(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                   IncomprFlowParam<T> const& parameters,
                   PhysUnits3D<T> const& units,
                   std::string const &fname, plint const iT)
{
  static bool init = true;

  plb_ofstream of;
  
  if(init){
    of.open(fname.c_str());
    of << "# time p_top p_bottom delta_p\n";
    init = false;
  } else{
    of.open(fname.c_str(),std::ostream::out|std::ostream::app);
  }

  T const time = units.getPhysTime(iT);

  pluint const nx = lattice.getNx(),
    ny = lattice.getNy(),
    nz = lattice.getNz();

  pluint const N = parameters.getResolution();
  Box3D const topPressureArea(N,nx-N,N,ny-N,nz-2*N,nz-2*N);
  Box3D const bottomPressureArea(N,nx-N,N,ny-N,N,N);

  T const rhoTop = computeAverage(*computeDensity(lattice,topPressureArea));
  T const rhoBottom = computeAverage(*computeDensity(lattice,bottomPressureArea));

  T const pTop = units.getPhysPress(rhoTop);
  T const pBottom = units.getPhysPress(rhoBottom);

  of << time << " " << pTop << " " << pBottom << " " << pBottom - pTop << std::endl;
  of.close();
}


int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);

    plint N(0);
    
    T u_in(0), u_max_re(0), uMax(0), dt_dem_target(1);
    std::string outDir;

    T cohesion_energy_density(0);
    
    InputParser input(argc,argv);

    T const uMaxReFactor = 5.;

    try{
      N = input.readInt("-N");
      u_in = input.readDouble("--u_in");
      u_max_re = uMaxReFactor*u_in;
      uMax = input.readDouble("--u_max_lb");
      outDir = input.readString("--outdir");
    } catch (std::exception &e) {
      pcout << "argument " << e.what() << " not given.\n";
      return -1;
    }

    if(input.cmdOptionExists("--dt_dem_max"))
      dt_dem_target = input.readDouble("--dt_dem_max");

    T rho_f(0.), rho_s(0.), nu_f(0.), r_(0.), volFrac(0.);
    
    try{
      rho_f = input.readDouble("--rho_f");
      rho_s = input.readDouble("--rho_s");
      nu_f = input.readDouble("--nu_f");
      r_ = input.readDouble("--radius");
      volFrac = input.readDouble("--vol_frac");
      if(input.cmdOptionExists("--cohesion_energy_density"))
        cohesion_energy_density = input.readDouble("--cohesion_energy_density");
    } catch(std::exception &e) {
      pcout << "material parameter error in " << e.what() << "\n";
      return -1;
    }

    T lx(0.), ly(0.), lz(0.), inlet_distance(0.);

    try{
      lx = input.readDouble("--lx");
      ly = input.readDouble("--ly");
      lz = input.readDouble("--lz");
      inlet_distance = input.readDouble("--inlet_distance");
    } catch(std::exception &e) {
      pcout << "geometry parameter error in " << e.what() << "\n";
      return -1;
    }

    T maxT(0.),writeT(0.),rampTime(0.);
    bool writeFlowfield(false);

    try{
      maxT = input.readDouble("--max_t");
      writeT = input.readDouble("--write_t");
      if(input.cmdOptionExists("--ramp_t"))
        rampTime = input.readDouble("--ramp_t");
      if(input.cmdOptionExists("--write_flowfield"))
        writeFlowfield = input.readBool("--write_flowfield");
    } catch(std::exception &e) {
      pcout << "time settings problem at " << e.what() << "\n";
      return -1;
    }

    std::string lbOutDir(outDir), demOutDir(outDir);
    lbOutDir.append("tmp/"); demOutDir.append("post/");
    global::directories().setOutputDir(lbOutDir);

    std::string argv_fname(outDir);
    argv_fname.append("argv.txt");
    plb_ofstream argv_file(argv_fname.c_str());
    for(int i=1;i<argc;i+=2)
      argv_file << argv[i] << " " << argv[i+1] << "\n";
    argv_file.close();

    
    LiggghtsCouplingWrapper wrapper(argv,global::mpi().getGlobalCommunicator());

    wrapper.setVariable("dmp_dir",demOutDir);

    wrapper.setVariable("lx",lx);
    wrapper.setVariable("ly",ly);
    wrapper.setVariable("lz",lz);
    wrapper.setVariable("d_inlet",inlet_distance);
    wrapper.setVariable("r_part",r_);
    wrapper.setVariable("v_frac",volFrac);

    wrapper.setVariable("e_cohesion",cohesion_energy_density);
    
    wrapper.execFile("in.lbdem");

    PhysUnits3D<T> units(2.*r_,u_max_re,nu_f,lx,ly,lz,N,uMax,rho_f);

    IncomprFlowParam<T> parameters(units.getLbParam());

    plint nx = parameters.getNx(), ny = parameters.getNy(), nz = parameters.getNz();

    // get lattice decomposition from LIGGGHTS and create lattice according to parallelization
    // given in the LIGGGHTS input script
    LatticeDecomposition lDec(parameters.getNx(),parameters.getNy(),parameters.getNz(),
                              wrapper.lmp);
    SparseBlockStructure3D blockStructure = lDec.getBlockDistribution();
    ExplicitThreadAttribution* threadAttribution = lDec.getThreadAttribution();
    plint envelopeWidth = 1;

    MultiBlockLattice3D<T, DESCRIPTOR> 
      lattice (MultiBlockManagement3D (blockStructure, threadAttribution, envelopeWidth ),
               defaultMultiBlockPolicy3D().getBlockCommunicator(),
               defaultMultiBlockPolicy3D().getCombinedStatistics(),
               defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
               new DYNAMICS );
    
    lattice.initialize();

    if( input.cmdOptionExists("--use_smagorinsky")
       && input.readBool("--use_smagorinsky") ){
      T cSmago = 0.15;
      if(input.cmdOptionExists("--c_smago"))
        cSmago = input.readDouble("--c_smago");
      
      defineDynamics(lattice,lattice.getBoundingBox(),new SMAGO_DYNAMICS);
    } else {
      defineDynamics(lattice,lattice.getBoundingBox(),new DYNAMICS);
    }

    lattice.periodicity().toggle(0,true);
    lattice.periodicity().toggle(1,true);
    lattice.periodicity().toggle(2,false);
      
    writeLogFile(parameters, "double periodic fluidized bed");
    
    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition =
      createLocalBoundaryCondition3D<T,DESCRIPTOR>();
    Box3D bottom(0,nx-1,0,ny-1,0,0), top(0,nx-1,0,ny-1,nz-1,nz-1);
    // Box3D bottom(1,nx-2,1,ny-2,0,0), top(1,nx-2,1,ny-2,nz-1,nz-1);
    // Box3D left(0,0,0,ny-1,0,nz-1),right(nx-1,nx-1,0,ny-1,0,nz-1);
    // Box3D front(1,nx-2,0,0,0,nz-1),back(1,nx-2,ny-1,ny-1,0,nz-1);

    // defineDynamics(lattice,left,new BounceBack<T,DESCRIPTOR>());
    // defineDynamics(lattice,right,new BounceBack<T,DESCRIPTOR>());
    // defineDynamics(lattice,front,new BounceBack<T,DESCRIPTOR>());
    // defineDynamics(lattice,back,new BounceBack<T,DESCRIPTOR>());

    boundaryCondition->setVelocityConditionOnBlockBoundaries(lattice,bottom);
    boundaryCondition->setPressureConditionOnBlockBoundaries(lattice,top);

    setBoundaryVelocity(lattice,bottom,Array<T,3>(0.,0.,0.));

    setBoundaryDensity(lattice,top,1.);

    const plint maxSteps = units.getLbSteps(maxT);
    const plint writeSteps = max<plint>(units.getLbSteps(writeT),1);

    T dt_phys = units.getPhysTime(1);
    plint demSubsteps = ceil(dt_phys/dt_dem_target);
    if(demSubsteps < 5) demSubsteps = 5;
    T dt_dem = dt_phys/(T)demSubsteps;
    plint dumpSteps = writeSteps*demSubsteps;

    pcout << "------------------------------\n"
          << "omega: " << parameters.getOmega() << "\n" 
          << "dt_phys: " << dt_phys << "\n"
	  << "dt_dem: " << dt_dem << "\n"
	  << "dt_dem_target: " << dt_dem_target << "\n"
	  << "substeps: " << demSubsteps << "\n"
          << "maxT: " << maxT << " | maxSteps: " << maxSteps << "\n"
          << "u_in: " << u_in << "\n"
      // Reynolds number requested by user is different from what is 
      // internally for unit conversion
          << "Re : " << parameters.getRe()/uMaxReFactor << "\n"
          << "writeT: " << writeT << " | writeSteps: " << writeSteps 
          << " | writeFlowfield " << writeFlowfield << "\n"
          << "grid size: " << nx << " " << ny << " " << nz << "\n"
          << "------------------------------" << std::endl;

    // set timestep and output directory
    wrapper.setVariable("t_step",dt_dem);
    wrapper.setVariable("dmp_stp",dumpSteps);

    T const gRed = (rho_s-rho_f)/rho_s*9.81;
    T const t_settle = sqrt(2*lz/gRed);
    plint const settleStp = ceil(2.*t_settle/1e-5);
    plint const settleStpReal = ceil(double(settleStp)/double(dumpSteps))*dumpSteps;
    wrapper.execCommand("timestep 1e-5");
    wrapper.execCommand("thermo 1000");
    wrapper.execCommand("thermo_style custom step atoms ke cpu");
    wrapper.execCommand("fix ts_check all check/timestep/gran 1000 0.01 0.01");
    wrapper.runUpto(settleStpReal-1+demSubsteps);
    
    wrapper.execFile("in2.lbdem");


    // preparations for ramping of velocity
    T current_uIn(0.), uIn_lb(uMax * u_in / u_max_re), rampInc(0.);
    plint rampSteps(0),nRamp(1000),singleRampStep(0);

    if(rampTime > 0.){
      rampSteps = units.getLbSteps(rampTime);
      singleRampStep = units.getLbSteps(rampTime/(T)nRamp);
      rampInc = uIn_lb / (T)nRamp;

      pcout << "rampSteps: " << rampSteps << "\n"
            << "singleRampStep: " << singleRampStep << std::endl;

    } else{
      current_uIn = uIn_lb;
      pcout << "no velocity ramping requested" << std::endl;
    }

    Array<plint,6> numSpongeCells;

    numSpongeCells[0] = 0;
    numSpongeCells[1] = 0;
    numSpongeCells[2] = 0;
    numSpongeCells[3] = 0;
    numSpongeCells[4] = 0;
    numSpongeCells[5] = 2*N;

    std::vector<MultiBlock3D*> args;
    args.push_back(&lattice);
  
    Box3D spongeBox(1,nx-2,1,ny-2,nz-1-numSpongeCells[5],nz-1);
    applyProcessingFunctional(new ViscositySpongeZone<T,DESCRIPTOR>(nx, ny, nz,
                                                                    parameters.getOmega(),
                                                                    numSpongeCells),
                              spongeBox,
                              args);
    

    clock_t start = clock();
    clock_t loop = clock();
    clock_t end = clock(); 
    // Loop over main time iteration.


    lattice.toggleInternalStatistics(false);

    std::string inlet_vel_fname(lbOutDir);
    inlet_vel_fname.append("u_in.txt");
    plb_ofstream invelfile(inlet_vel_fname.c_str());
    invelfile << "# t u_in" << std::endl;
    invelfile.close();

    std::string pressure_fname(lbOutDir);
    pressure_fname.append("pressure.txt");

    for (plint iT=0; iT<=maxSteps; ++iT) {
      
      if(singleRampStep > 0 && iT%singleRampStep == 0){
        if(current_uIn < uIn_lb){
          current_uIn += rampInc;
          if(current_uIn > uIn_lb) current_uIn = uIn_lb;
          setBoundaryVelocity(lattice,bottom,Array<T,3>(0.,0.,current_uIn));
          pcout << "increased inlet velocity to " << units.getPhysVel(current_uIn) << std::endl;
          plb_ofstream invelfile(inlet_vel_fname.c_str(),std::ostream::out|std::ostream::app);
          invelfile << units.getPhysTime(iT) << " " << units.getPhysVel(current_uIn) << std::endl;
          invelfile.close();
        }
      }

      pcout << "setSpheres" << std::endl;
      bool initWithVel = false;
      setSpheresOnLattice(lattice,wrapper,units,initWithVel);
      

      if(iT%writeSteps == 0 && iT > 0){ // LIGGGHTS does not write at timestep 0
        if(writeFlowfield)
          writeVTK(lattice,parameters,units,iT);
        writePressure(lattice,parameters,units,pressure_fname,iT);
      }
      
      pcout << "collideAndStream" << std::endl;
      lattice.collideAndStream();

      // T const uMax_current = computeMax(*computeVelocityNorm(lattice,lattice.getBoundingBox()));

      // pcout << "maximum velocity in domain: " << uMax_current << std::endl;

      pcout << "getForcesFromLattice" << std::endl;
      getForcesFromLattice(lattice,wrapper,units);

      pcout << "runDEM" << std::endl;
      wrapper.run(demSubsteps);
      

      end = clock();
      T time = difftime(end,loop)/((T)CLOCKS_PER_SEC);
      T totaltime = difftime(end,start)/((T)CLOCKS_PER_SEC);
      T mlups = ((T) (lattice.getNx()*lattice.getNy()*lattice.getNz()))/time/1e6;
      pcout << "time: " << time << " " ;
      pcout << "calculating at " << mlups << " MLU/s"
            << " | total time running: " << totaltime << std::endl;
      loop = clock();
    
    }
    T totaltime = difftime(end,start)/((T)CLOCKS_PER_SEC);
    T totalmlups = ((T) (lattice.getNx()*lattice.getNy()*lattice.getNz()*(maxSteps+1)))/totaltime/1e6;
    pcout << " ********************** \n"
          << "total time: " << totaltime
          << " calculating at " << totalmlups << " MLU/s" << std::endl;

}
