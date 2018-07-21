
///////////////////////////////////////////////////////////////////////////////
//
// File UnsteadyAdvectionDiffusion.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Unsteady advection-diffusion solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <ADRSolver/EquationSystems/UnsteadyAdvectionDiffusion.h>

using namespace std;

namespace Nektar
{
    
    double ADPtcount = 0, gammaeta;
    Array<OneD, NekDouble> coord_0, coord_1, coord_2, 
        interpCoord_x, interpCoord_y, interpCoord_z, PboldInterp, 
        kadh, kadhNonScaled, Bkadh, rtemp, cover, tfAvail, 
        ze9btot, e2btot, ze5btot, ze10btot, ze8btot, z2btot, 
        ze9btotInterp, e2btotInterp, ze5btotInterp, ze10btotInterp, ze8btotInterp, 
        z2btotInterp, 
        PLPba, PLPsea;

    //sigmar holds all future information
    Array< OneD, Array<OneD, NekDouble> >  sigmar, outArrayReact, fpReactInterp, fp, fpInterp, pb, se, se2, pb2;

    Array<OneD, NekDouble> geta;       
    //sigmaz is the CURRENT release value
    Array<OneD, NekDouble>  sigmaz;
    Array<OneD, NekDouble> Pbnew, Pbold, NonDpb, NonDse, phiT;
    std::map<std::string, int> varFPNames, varSENames, varPBNames;

    // Reaction zone
    double m_rzl, m_rzr, m_rzt, nonDAdvVel; 
    // Bottom row and reaction zone andleft boundary coordinates  
    Array<OneD, int> brCoord, rzCoord, ofCoord;

    int nScaledDomainPts = 0, nSolutionPts = 0; //no. of points in scaled domain and domain
    int nScaledBtr = 0; //no. of points along bot row scaled
    int nBtr = 0; //no. of points along bottom row scaled 
    int numRZ = 0; //no. of points along rz scaled
    int nOutflow = 0, nScOutflow = 0; //no. of points along outflow boundary and scaled outflow boundary

    NekDouble prev = 0, check = 0; //output PB and SE variables to txt files
    NekDouble OneDptscale; //scaling factor for scaled domain (fp and Platelets)
    NekDouble OneDptscaleS; //scaling factor for scaled domain (pb and se)

    string UnsteadyAdvectionDiffusion::className
    = SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "UnsteadyAdvectionDiffusion",
        UnsteadyAdvectionDiffusion::create);
    
    UnsteadyAdvectionDiffusion::UnsteadyAdvectionDiffusion(
        const LibUtilities::SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pSession),
          AdvectionSystem(pSession)
    {
        m_planeNumber = 0;
    }
    
    /**
     * @brief Initialisation object for the unsteady linear advection 
     * diffusion equation.
     */
    void UnsteadyAdvectionDiffusion::v_InitObject()
    {
        AdvectionSystem::v_InitObject();
	
        initPAC();

        m_session->LoadParameter("wavefreq",   m_waveFreq, 0.0);

        m_session->LoadParameter("epsilon",    m_epsilon,  0.0);
        
        // turn on substepping
        m_session->MatchSolverInfo("Extrapolation", "SubStepping",
                                   m_subSteppingScheme, false);
        
        // Define Velocity fields
        m_velocity = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
        std::vector<std::string> vel;
        vel.push_back("Vx");
        vel.push_back("Vy");
        vel.push_back("Vz");
        vel.resize(m_spacedim);
        
        EvaluateFunction(vel, m_velocity, "AdvectionVelocity");
        
        m_session->MatchSolverInfo(
            "SpectralVanishingViscosity", "True", m_useSpecVanVisc, false);
        
        if(m_useSpecVanVisc)
        {
            m_session->LoadParameter("SVVCutoffRatio",m_sVVCutoffRatio,0.75);
            m_session->LoadParameter("SVVDiffCoeff",m_sVVDiffCoeff,0.1);
        }        

        // Type of advection and diffusion classes to be used
        switch(m_projectionType)
        {
            // Discontinuous field 
        case MultiRegions::eDiscontinuous:
        {
            // Do not forwards transform initial condition
            m_homoInitialFwd = false;

            // Advection term
            string advName;
            string riemName; 
            m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
            m_advObject = SolverUtils::GetAdvectionFactory().
                CreateInstance(advName, advName);
            m_advObject->SetFluxVector(&UnsteadyAdvectionDiffusion::
                                       GetFluxVectorAdv, this);
            m_session->LoadSolverInfo("UpwindType", riemName, "Upwind");
            m_riemannSolver = SolverUtils::GetRiemannSolverFactory().
                CreateInstance(riemName);
            m_riemannSolver->SetScalar("Vn", &UnsteadyAdvectionDiffusion::
                                       GetNormalVelocity, this);
            m_advObject->SetRiemannSolver(m_riemannSolver);
            m_advObject->InitObject      (m_session, m_fields);
                
            // Diffusion term
            std::string diffName;
            m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
            m_diffusion = SolverUtils::GetDiffusionFactory().
                CreateInstance(diffName, diffName);
            m_diffusion->SetFluxVector(&UnsteadyAdvectionDiffusion::
                                       GetFluxVectorDiff, this);
            m_diffusion->InitObject(m_session, m_fields);

            ASSERTL0(m_subSteppingScheme == false,"SubSteppingScheme is not set up for DG projection");
            break;
        }
        // Continuous field 
        case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
        {
            // Advection term
            std::string advName;
            m_session->LoadSolverInfo("AdvectionType", advName, 
                                      "NonConservative");
            m_advObject = SolverUtils::GetAdvectionFactory().
                CreateInstance(advName, advName);
            m_advObject->SetFluxVector(&UnsteadyAdvectionDiffusion::
                                       GetFluxVectorAdv, this);

            if(advName.compare("WeakDG") == 0)
            {
                string riemName;
                m_session->LoadSolverInfo("UpwindType", riemName, "Upwind");
                m_riemannSolver = SolverUtils::GetRiemannSolverFactory().
                    CreateInstance(riemName);
                m_riemannSolver->SetScalar("Vn",
                                           &UnsteadyAdvectionDiffusion::
                                           GetNormalVelocity, this);
                m_advObject->SetRiemannSolver(m_riemannSolver);
                m_advObject->InitObject      (m_session, m_fields);
            }

            // In case of Galerkin explicit diffusion gives an error
            if (m_explicitDiffusion)
            {
                ASSERTL0(false, "Explicit Galerkin diffusion not set up.");
            }
            // In case of Galerkin implicit diffusion: do nothing
            break;
        }
        default:
        {
            ASSERTL0(false, "Unsupported projection type.");
            break;
        }
        }

        m_ode.DefineImplicitSolve (
            &UnsteadyAdvectionDiffusion::DoImplicitSolve, this);

        if(m_subSteppingScheme) // Substepping
        {
            ASSERTL0(m_projectionType == MultiRegions::eMixed_CG_Discontinuous,
                     "Projection must be set to Mixed_CG_Discontinuous for "
                     "substepping");
            SetUpSubSteppingTimeIntegration(
                m_intScheme->GetIntegrationMethod(), m_intScheme);

        }
        else // Standard velocity correction scheme
        {
            m_ode.DefineOdeRhs(&UnsteadyAdvectionDiffusion::DoOdeRhs, this);
        }

        if (m_projectionType == MultiRegions::eDiscontinuous &&
            m_explicitDiffusion == 1)
        {
            m_ode.DefineProjection(&UnsteadyAdvectionDiffusion::DoOdeProjection, this);
        }
    }

    /**
     * @ brief Platelet Aggregation and Coagulation model (PAC) initializer
     */
    void UnsteadyAdvectionDiffusion::initPAC()
    {
        // Number of fields (variables of the problem)
        int nVariables = m_fields.num_elements();
        
        // Number of solution points
        nSolutionPts = GetNpoints();

        OneDptscale = 1.5;
        OneDptscaleS = 1.0;

        gammaeta = m_epsilon*10/pow(3e-4, 2);
        //cout<<"\n gammatea init to "<<gammaeta<<" m_epsilon ="<<m_epsilon;
        gammaeta = 55.5555555556;
        //cout<<"\n gammaeta modified to "<<gammaeta;
        // Initialisation of higher-space variables      
        nScaledDomainPts = m_fields[0]->Get1DScaledTotPoints(OneDptscale);
        //cout<<"\nnScaledDomainPts = "<<nScaledDomainPts<<" and nSolutionPts="<<nSolutionPts<<endl;

        nonDAdvVel = 0.1; //ND-AdvCoefficient

        pb = Array<OneD, Array<OneD, NekDouble> >(numPBVars);
        se = Array<OneD, Array<OneD, NekDouble> >(numSEVars);

        geta = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);

        // Initialize total variables
        ze9btotInterp =  Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
        ze9btot =  Array<OneD, NekDouble>(nSolutionPts, 0.0);
        e2btotInterp =  Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
        e2btot =  Array<OneD, NekDouble>(nSolutionPts, 0.0);
        ze5btotInterp =  Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
        ze5btot =  Array<OneD, NekDouble>(nSolutionPts, 0.0);
        ze10btotInterp =  Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
        ze10btot =  Array<OneD, NekDouble>(nSolutionPts, 0.0);
        ze8btotInterp =  Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
        ze8btot =  Array<OneD, NekDouble>(nSolutionPts, 0.0);
        z2btotInterp =  Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
        z2btot =  Array<OneD, NekDouble>(nSolutionPts, 0.0);
        
        PLPba =  Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
        PLPsea =  Array<OneD, NekDouble>(nScaledDomainPts, 0.0);

        // Initialize the vectors for PB variables
        for (int i = 0; i < numPBVars; i++ )
            pb[i] = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);

        // Initialize the vectors for SE variables
        for (int i = 0; i < numSEVars; i++ )
            se[i] = Array<OneD, NekDouble>(nScaledDomainPts, 0.0); 

        
        rtemp = Array<OneD, NekDouble>(rlen, 0.0);
        sigmar = Array<OneD, Array<OneD, NekDouble> >(nScaledDomainPts);
        for( int i = 0; i < nScaledDomainPts; i++)
            sigmar[i] = Array<OneD, NekDouble>(rlen,0.0);
        sigmaz = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
        Pbnew = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
        Pbold = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);

        // initialize non-dimentionalization parameters 
        NonDse = Array<OneD, NekDouble>(numSEVars);
        NonDpb = Array<OneD, NekDouble>(numPBVars); 


        // setup rtemp for ADP release info
        double ir = 0.0;
        for(int i = 0; i<rlen; i++)
        {
            double tsigma = (i - 1)*5.0/(rlen - 1);
            rtemp[i] = exp( -pow((tsigma-3.0),2) );
            ir = ir + rtemp[i]*5.0/(rlen - 1);
        }
        Vmath::Smul(rlen, 1/ir, rtemp, 1, rtemp, 1);
	
        initMaps();
        
        fp = Array<OneD, Array<OneD, NekDouble> >(nVariables);
        fpInterp = Array<OneD, Array<OneD, NekDouble> >(nVariables);
        se2 = Array<OneD, Array<OneD, NekDouble> >(numSEVars);
        pb2 = Array<OneD, Array<OneD, NekDouble> >(numPBVars);

        for( int i = 0; i < nVariables; i++)
        {
            fp[i] = Array<OneD, NekDouble>(nSolutionPts, 0.0);
            fpInterp[i] = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
        }
        for( int i = 0; i < numPBVars; i++)
        {
            pb2[i] = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
        }
        for( int i = 0; i < numSEVars; i++)
        {
            se2[i] = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
        }

        PboldInterp = Array<OneD, NekDouble> (nScaledDomainPts);
        kadh = Array<OneD, NekDouble>(nScaledDomainPts, 0.0); //for RZ = kadhmax, else 0   
        kadhNonScaled = Array<OneD, NekDouble>(nSolutionPts, 0.0); //for RZ = kadhmax, else 0   
        cover = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);

        // get coefficients to determine the RZ
        coord_0 =  Array<OneD, NekDouble> (nSolutionPts);
        coord_1 =  Array<OneD, NekDouble> (nSolutionPts);
        coord_2 =  Array<OneD, NekDouble> (nSolutionPts); 
        interpCoord_x =  Array<OneD, NekDouble> (nScaledDomainPts);
        interpCoord_y =  Array<OneD, NekDouble> (nScaledDomainPts);
        interpCoord_z =  Array<OneD, NekDouble> (nScaledDomainPts);
        m_fields[0]->GetCoords(coord_0, coord_1, coord_2);
        m_fields[0]->PhysInterp1DScaled(OneDptscale, coord_0, interpCoord_x);
        m_fields[0]->PhysInterp1DScaled(OneDptscale, coord_1, interpCoord_y);
        m_fields[0]->PhysInterp1DScaled(OneDptscale, coord_2, interpCoord_z);

        // Interpolation to higher space of fields and initialization of interp_react
        outArrayReact = Array< OneD, Array<OneD, NekDouble> >(nVariables);
        fpReactInterp =  Array<OneD, Array<OneD, NekDouble> >(nVariables);

        for (int i = 0; i < nVariables; ++i)
        {
            fpInterp[i] = Array<OneD, NekDouble>(nScaledDomainPts);
              
            // also initialize the reactionInterp array:
            fpReactInterp[i] = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
            outArrayReact[i] = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
        }


        for (int i = 0; i < nSolutionPts; i++)
        {
            if(coord_1[i] == 0.0 )//y
                nBtr++;
            if(coord_0[i] == 1.0)//x
                nOutflow++;
         }   
        
        cout<<"\nBottom most row has "<<nBtr<<" points.\n";
        
        
        for (int i = 0; i < nScaledDomainPts; i++)
        {
            if(interpCoord_y[i] == 0.0 )//y
            {
                nScaledBtr++;
            }
            if(interpCoord_x[i] >= 0.999999999999999 )//x
            {
                nScOutflow++;
            }
        }   
        cout<<"\nBottom most row has interpolated"<<nScaledBtr<<" points.\n";
        cout<<"\nOutflow boundary has interpolated"<<nScOutflow<<" points.";
        brCoord =  Array<OneD, int> (nScaledBtr); //bot row index
        ofCoord = Array<OneD, int>(nScOutflow); //outflow col index
        
        for (int i = 0, ct1 = 0, ct2 = 0; i < nScaledDomainPts; i++)
        {
            if(interpCoord_y[i] == 0.0 )//y
            {
                brCoord[ct1++] = i;
            }
            if(interpCoord_x[i] >= 0.999999999999999 )//x
            {
                ofCoord[ct2++] = i;
            }
        }   
        // for (int i = 0, ct1 = 0; i<nSolutionPts; i++)
        //     if(coord_1[i] == 0.0)//y
        //         brCoord[ct1++] = i;
       
        // get reaction zone parameters 
        if (m_session->DefinesParameter("rzl"))
        {
            m_rzl = m_session->GetParameter("rzl");
        }
        if (m_session->DefinesParameter("rzr"))
        {
            m_rzr = m_session->GetParameter("rzr");
        }
        if (m_session->DefinesParameter("rzt"))
        {
            m_rzt = m_session->GetParameter("rzt");
        }

        // 2D: Loop over quadrature points to find points inside RZ 
        // and assign them kadh = kadhmax 
        for (int i = 0; i < nScaledDomainPts; i++)
        {
            if(interpCoord_x[i] >= m_rzl && interpCoord_x[i] <= m_rzr)
            {  
                if(interpCoord_y[i] == 0) numRZ++;
                if(interpCoord_y[i] >= 0 && interpCoord_y[i] <= m_rzt)
                {
                    kadh[i] = kadhmax;
                }
                
            }
        }
  
        // for(int i = 0; i<nSolutionPts; i++)
        // {
        //     if(coord_0[i]>=m_rzl && coord_0[i]<=m_rzr && coord_1[i]==0)
        //         numRZ++;
        // }
      
        rzCoord = Array<OneD, int>(numRZ);

        for (int i = 0, ct = 0; i < nScaledDomainPts; i++)
        {
            if(interpCoord_x[i] >= m_rzl && interpCoord_x[i] <= m_rzr && interpCoord_y[i] == 0)
//            if(coord_0[i]>=m_rzl && coord_0[i]<=m_rzr && coord_1[i]==0)

            {
                rzCoord[ct++] = i; //interpCoord_x[i];
            }
        }

        Bkadh = Array<OneD,NekDouble>(nScaledBtr);
        // Galerkin Project kadh back to non-scaled domain
        m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, kadh, kadhNonScaled);
        //check if kadh or kadhNonScaled is negative because of this projection
        for (int i = 0; i<nScaledDomainPts; i++)
            if(kadh[i]<0)
            { 
                cout<<"\n kadh["<<
                    i<<"]="<<kadh[i];
                exit(0);
            }
        cout<<"\n Scaled kadh is all non-negative!\
              Now checking projection effects on positivity";
        for (int i = 0; i<nSolutionPts; i++)
            if(kadhNonScaled[i]<0)
            { 
                cout<<"\n kadhNonScaled["<<
                    i<<"]="<<kadhNonScaled[i];
                exit(0);
            }

        Vmath::Gathr(nScaledBtr, kadhNonScaled, brCoord, Bkadh);

        // Initialize tfAvail vector
        tfAvail = Array<OneD, NekDouble>(nScaledBtr, 0.0);

        //se[SEtf] is initialized to 1 at rz boundary
        for(int i = 0; i<numRZ; i++)
            se[SEtf][rzCoord[i]] = 1; //Can use Gathr/Scatr


        prev = -1; // Output PB and SE variables in txt files
        //cout<<"\n Calling OutputPBandSE!";
        OutputPBandSE(0);
        check++;

    }

    /**
     * @ brief All variable name maps initializer
     */
    void UnsteadyAdvectionDiffusion::initMaps()
    {
        // Number of fields (variables of the problem)
        int nVariables = m_fields.num_elements();

        for (int i = 0; i<nVariables; i++)
        {
            varFPNames.insert( std::pair<std::string, int>(m_session->GetVariable(i), i));
            
        }

        varPBNames.insert( std::pair<std::string, int>("PBz5", 0));
        varPBNames.insert( std::pair<std::string, int>("PBz8e2", 1));
        varPBNames.insert( std::pair<std::string, int>("PBz5e2", 2));
        varPBNames.insert( std::pair<std::string, int>("PBe2", 3));
        varPBNames.insert( std::pair<std::string, int>("PBz5e10", 4));
        varPBNames.insert( std::pair<std::string, int>("PBe5", 5));
        varPBNames.insert( std::pair<std::string, int>("PBapce5", 6));
        varPBNames.insert( std::pair<std::string, int>("PBpro", 7));
        varPBNames.insert( std::pair<std::string, int>("PBz2pro", 8));
        varPBNames.insert( std::pair<std::string, int>("PBz10tenstar", 9));
        varPBNames.insert( std::pair<std::string, int>("PBz8e10", 10));
        varPBNames.insert( std::pair<std::string, int>("PBz10ten", 11));
        varPBNames.insert( std::pair<std::string, int>("PBe10", 12));
        varPBNames.insert( std::pair<std::string, int>( "PBz10", 13));
        varPBNames.insert( std::pair<std::string, int>("PBe8", 14));
        varPBNames.insert( std::pair<std::string, int>("PBtenstar", 15));
        varPBNames.insert( std::pair<std::string, int>("PBten", 16));
        varPBNames.insert( std::pair<std::string, int>("PBz8", 17));
        varPBNames.insert( std::pair<std::string, int>("PBapce8", 18));
        varPBNames.insert( std::pair<std::string, int>("PBz2", 19));
        varPBNames.insert( std::pair<std::string, int>("PBe9", 20));
        varPBNames.insert( std::pair<std::string, int>("PBe9star", 21));
        varPBNames.insert( std::pair<std::string, int>("PBz9", 22));

        varSENames.insert( std::pair<std::string, int>("SEz10e7", 0));
        varSENames.insert( std::pair<std::string, int>("SEe7", 1));
        varSENames.insert( std::pair<std::string, int>("SEtfpie10e7", 2));
        varSENames.insert( std::pair<std::string, int>("SEz7", 3));
        varSENames.insert( std::pair<std::string, int>("SEz7e10", 4));
        varSENames.insert( std::pair<std::string, int>("SEtf" , 5));
        varSENames.insert( std::pair<std::string, int>("SEz7e2", 6));
        varSENames.insert( std::pair<std::string, int>("SEz9e7", 7));
         
        NonDpb[PBz5] = s5*1.0e9;
        NonDpb[PBz8e2] = s8*1.0e9;
        NonDpb[PBz5e2] = s2s*1.0e9; 
        NonDpb[PBe2] = s2s*1.0e9;
        NonDpb[PBz5e10] = s10*1.0e9;
        NonDpb[PBe5] = s5*1.0e9;
        NonDpb[PBapce5] = s5*1.0e9;
        NonDpb[PBpro] = s5*1.0e9;
        NonDpb[PBz2pro] = s5*1.0e9;
        NonDpb[PBz10tenstar] = s8*1.0e9;
        NonDpb[PBz8e10] = s8*1.0e9;
        NonDpb[PBz10ten] = s8*1.0e9; 
        NonDpb[PBe10] = s10*1.0e9;
        NonDpb[PBz10] = s10*1.0e9;
        NonDpb[PBe8] = s8*1.0e9;
        NonDpb[PBtenstar] = s8*1.0e9;
        NonDpb[PBten] = s8*1.0e9;
        NonDpb[PBz8] = s8*1.0e9; 
        NonDpb[PBapce8] = s8*1.0e9;
        NonDpb[PBz2] = s2*1.0e9;
        NonDpb[PBe9] = s9*1.0e9;
        NonDpb[PBe9star] = s91*1.0e9;
        NonDpb[PBz9] = s9*1.0e9;
        
        NonDse[SEz10e7] = s7d*1e12 ;
        NonDse[SEe7] = s7d*1e12; 
        NonDse[SEtfpie10e7] = s7d*1e12;
        NonDse[SEz7] = s7d*1e12;
        NonDse[SEz7e10] = s7d*1e12;
        NonDse[SEtf] = 1;
        NonDse[SEz7e2] = s7d*1e12;
        NonDse[SEz9e7] = s7d*1e12;

    }

    /**
     * @brief Unsteady linear advection diffusion equation destructor.
     */
    UnsteadyAdvectionDiffusion::~UnsteadyAdvectionDiffusion()
    {
    }
    
    /**
     * @brief Get the normal velocity for the unsteady linear advection 
     * diffusion equation.
     */
    Array<OneD, NekDouble> &UnsteadyAdvectionDiffusion::GetNormalVelocity()
    {
        return GetNormalVel(m_velocity);
    }


    Array<OneD, NekDouble> &UnsteadyAdvectionDiffusion::GetNormalVel(
        const Array<OneD, const Array<OneD, NekDouble> > &velfield)
    {
        // Number of trace (interface) points
        int i;
        int nTracePts = GetTraceNpoints();

        // Auxiliary variable to compute the normal velocity
        Array<OneD, NekDouble> tmp(nTracePts);
        m_traceVn = Array<OneD, NekDouble>(nTracePts, 0.0);

        // Reset the normal velocity
        Vmath::Zero(nTracePts, m_traceVn, 1);

        for (i = 0; i < velfield.num_elements(); ++i)
        {
            m_fields[0]->ExtractTracePhys(velfield[i], tmp);

            Vmath::Vvtvp(nTracePts,
                         m_traceNormals[i], 1,
                         tmp, 1,
                         m_traceVn, 1,
                         m_traceVn, 1);
        }
        
        return m_traceVn;
    }

    NekDouble NekTanh(NekDouble d) {return std::tanh(d);}

    
    /**
     * @brief Compute the right-hand side for the unsteady linear advection 
     * diffusion problem.
     * 
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void UnsteadyAdvectionDiffusion::DoOdeRhs(
        const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
        Array<OneD,        Array<OneD, NekDouble> >&outarray,
        const NekDouble time)
    {

        //cout<<"\n******DOODERHS*************";
        // Number of fields (variables of the problem)
        int nVariables = inarray.num_elements();
        
        // Number of solution points
        int nSolutionPts = GetNpoints();
      
    	// Number of dimensions
    	int nDim = m_velocity.num_elements(); 
	
    	//edited by vidhi
    
        // for(int i = 0; i<inarray.num_elements(); i++)
        // for(int k = 0; k<nSolutionPts; k++)        
        // if(inarray[i][k]<0) {cout<<"\n Advection inarray got negative val for i = "<<i; exit(0);};


        double xchar, tchar;
	
        xchar = m_session->GetParameter("xchar");
        tchar = m_session->GetParameter("tchar");

        // Update Pbnew and Pbold
	
        Array<OneD, NekDouble>  Ptemp(nScaledDomainPts);  

        Vmath::Vadd(nScaledDomainPts, PLPba, 1, PLPsea, 1, Pbnew, 1);
        Vmath::Vsub(nScaledDomainPts, Pbnew, 1, Pbold, 1, Pbnew, 1);
        Vmath::Vadd(nScaledDomainPts, PLPba, 1, PLPsea, 1, Pbold, 1);	
	
        ADPtcount = ADPtcount + tchar; // dt = tchar ? 
        //if(ADPtcount > 0.25/tchar) //UPDATE EACH TIME-STEP
        //{
        //cout<<"Updating ADP release at time "<<time<<"\n";
        UpdateADPRelease();
        ADPtcount = 0;
        Pbnew = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
        //}
        // Find WphiT for Pmu and Pma:
        
        phiT = Array<OneD, NekDouble>(nSolutionPts);  

        /*testing..remove after verification
        ofstream myfile;
        std::string filename =  "PmuInit.csv.0";
        myfile.open(filename.c_str());
        myfile << "x coord;y coord;z coord;scalar\n";
        for(int i = 0; i < nSolutionPts; i++)
            myfile<<coord_0[i]<<";"<<coord_1[i]<<";"<<coord_2[i]<<";"<<inarray[Pmu][i]<<"\n";//
            myfile.close(); exit(0); verified left border = 1*/
        ////////

        Vmath::Vadd(nSolutionPts, inarray[Pmu], 1, 
                    inarray[Pma], 1, phiT, 1);
        Vmath::Smul(nSolutionPts, phitfac, 
                    phiT, 1, phiT, 1);
        //cout<<"\n max(pmu*Pma*phitfac) = "<<Vmath::Vmax(nSolutionPts, phiT, 1);
        Vmath::Vadd(nSolutionPts, PLPsea, 1, 
                    phiT, 1, phiT, 1);
        Vmath::Vadd(nSolutionPts, PLPba, 1, 
                    phiT, 1, phiT, 1);
        if(Vmath::Vmax(nSolutionPts, phiT,1) > 1 ) 
            cout<<"\tPlatelet sum exceeds volume max !!! Error!!! value="<<
                Vmath::Vmax(nSolutionPts, &phiT[0],1)<<"\t phitfac = "<<phitfac;
        //cout<<"\n Max phiT="<<Vmath::Vmax(nSolutionPts, &phiT[0],1);
        Vmath::Neg(nSolutionPts, phiT, 1);
        Vmath::Sadd(nSolutionPts, 1.0, phiT, 1, 
                    phiT, 1.0); /* 1-phiT*/
        Vmath::Smul(nSolutionPts, PI, 
                    phiT, 1, phiT, 1); /*PI*(1-phiT)*/
	      
        transform(phiT.begin(), phiT.end(), phiT.begin(), NekTanh);/*WphiT*/
	
        Array<OneD, Array<OneD, NekDouble> > newAdvVel(nDim);
        Array<OneD, Array<OneD, NekDouble> > mobileAdvVel(nDim);
    	for (int i = 0; i < nDim; i++)
    	{
            if(i == 0)
                newAdvVel[i] = Array<OneD, NekDouble>(nSolutionPts, nonDAdvVel);
            else
                newAdvVel[i] = Array<OneD, NekDouble>(nSolutionPts, 0.0);
            mobileAdvVel[i] = Array<OneD, NekDouble>(nSolutionPts, 0.0);
        }

        // Multiply newAdvVel with WphiT
        Vmath::Vmul(nSolutionPts, newAdvVel[0], 1, phiT, 1, mobileAdvVel[0], 1);

        // if velocity in y-direction non-zero, multiply newAdvvel[1] to WphiT as well
        for(int i = 0; i<nSolutionPts; i++)
          if(mobileAdvVel[0][i] < 0.0) 
          {
              mobileAdvVel[0][i] = 0.0;
              cout<<"\n mobileAdvvel[0]["<<i<<"] is negative;";
              exit(0);
          }
        
        Array<OneD, boost::shared_ptr<Nektar::MultiRegions::ExpList> > newm_fields(2); //only Pmu and Pma
        newm_fields[Pmu] = m_fields[Pmu];
        newm_fields[Pma] = m_fields[Pma];
        Array<OneD, Array<OneD, NekDouble> > newinarray(2); //only Pmu and Pma
        newinarray[Pmu] = Array<OneD, NekDouble>(nSolutionPts);
        newinarray[Pma] = Array<OneD, NekDouble>(nSolutionPts);
        //newinarray[0] = inarray[Pmu];
        Vmath::Vcopy(nSolutionPts, inarray[Pmu], 1, newinarray[Pmu], 1);
        
        //newinarray[1] = inarray[Pma];
        Vmath::Vcopy(nSolutionPts, inarray[Pma], 1,  newinarray[Pma], 1);
        Array<OneD, Array<OneD, NekDouble> > newoutarray(2); //only Pmu and Pma
        newoutarray[0] = Array<OneD, NekDouble>(nSolutionPts);//outarray[Pmu];
        newoutarray[1] = Array<OneD, NekDouble>(nSolutionPts);//outarray[Pma];
        
        /*cout<<"\nBefore Advection. \nPmuMax= "<<Vmath::Vmax(nSolutionPts, newinarray[Pmu],1); 
          cout<<" Min="<<Vmath::Vmin(nSolutionPts, newinarray[Pmu],1) <<"\t";          
          cout<<"\nChecking mobileAdvvel: Max="<<Vmath::Vmax(nSolutionPts, &mobileAdvVel[0][0],1)
          <<" min = "<<Vmath::Vmin(nSolutionPts, &mobileAdvVel[0][0],1)<<endl;*/
        
        
        // RHS computation using the new advection base class (only for Pmu and Pma)
        m_advObject->Advect(2,  newm_fields, mobileAdvVel, newinarray, newoutarray, time);
        Vmath::Vcopy(nSolutionPts, newoutarray[0], 1, outarray[Pmu],1);
        Vmath::Vcopy(nSolutionPts, newoutarray[1], 1, outarray[Pma],1);
        
        
        for(int i = 0; i<nSolutionPts; i++)
        {
            if(outarray[Pmu][i] < 0.0) outarray[Pmu][i] = 0.0;
            if(outarray[Pma][i] < 0.0) outarray[Pma][i] = 0.0;
        }
        /*cout<<"time="<<time;*/
        
        
        // only ADP through nVariables
        Array<OneD, boost::shared_ptr<Nektar::MultiRegions::ExpList> > newm_fields2(nVariables-ADP);
        Array<OneD, Array<OneD, NekDouble> > newinarray2(nVariables-ADP);    
        Array<OneD, Array<OneD, NekDouble> > newoutarray2(nVariables-ADP);
        for(int i = ADP, ct = 0; i<nVariables; i++, ct++)
        {
            newm_fields2[ct] = m_fields[i]; 
            newinarray2[ct] = inarray[i];
            newoutarray2[ct] = outarray[i];
        }
        
        // RHS computation using the new advection base class (only for ADP through nVariables)
        m_advObject->Advect(nVariables-ADP,  newm_fields2, newAdvVel, newinarray2, newoutarray2, time);
        for(int i = ADP, ct = 0; i<nVariables; i++, ct++)
        {
            Vmath::Vcopy(nSolutionPts, newoutarray2[ct], 1, outarray[i],1);
        }

        //cout<<"\nDone Advection. ";

        //Negate the RHS (only the ones that advect)
        /*for( int i = 0; i<nVariables; i++)
        {
            if(i !=  Peta)
                Vmath::Neg(nSolutionPts, outarray[i], 1); 
        }*/
	   
        //check if outarray is +ve:
        for( int i = 0; i<nVariables; i++)
        {
            for(int j = 0; j<nSolutionPts; j++)
            {
                if(abs(outarray[i][j])<1e-14)
                    outarray[i][j]=0;
                if(outarray[i][j]<0) 
                {
                    cout<<"\n outarray was negative after advection for i = "<<
                        i<<" and j = "<<
                        j<<" val = "<<
                        outarray[i][j];
                }
            }
         }
        

        // Reaction: All Platelets and chemicals
         
        DoReact(inarray, outarray, time);

        /*
        //check if outarray is +ve:
        for( int i = 0; i<nVariables; i++)
        {
            for(int j = 0; j<nSolutionPts; j++)
                if(outarray[i][j]<0) 
                {
                    cout<<"\n outarray was negative after DoReact for i = "<<
                        i<<" and j = "<<
                        j<<" val = "<<
                        outarray[i][j];
                }
         }*/
       
        // No explicit diffusion for CG
        if (m_projectionType == MultiRegions::eDiscontinuous)
        {
            m_diffusion->Diffuse(nVariables, m_fields, inarray, outarray);

            for (int i = 0; i < nVariables; ++i)
            {
                Vmath::Vadd(nSolutionPts, &outarray[i][0], 1, 
                            &outarray[i][0], 1, &outarray[i][0], 1);
            }
        }
        
    }

/**
 * @brief Compute the right-hand side reaction part for the unsteady linear advection 
 * diffusion reaction problem.
 * 
 * @param inarray    Given fields.
 * @param outarray   Calculated solution.
 * @param outArrayReact   Calculated reactions.
 */
    void UnsteadyAdvectionDiffusion::DoReact(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        Array<OneD, Array<OneD, NekDouble> > &outarray,
        NekDouble time)
    {
        // Number of fields (variables of the problem)
        int nVariables = inarray.num_elements();

        //cout<<"\n******DOREACT************* at "<<time;

        // Number of solution points
        int nSolutionPts = GetNpoints();   

        // characteristic time (dt)
        double tchar = m_session->GetParameter("tchar"); 
        
        //for RK-2
        Array<OneD, Array<OneD, NekDouble> >k1 = Array<OneD, Array<OneD, NekDouble> >(numSEVars+numPBVars+2);
        Array<OneD, Array<OneD, NekDouble> >k2 = Array<OneD, Array<OneD, NekDouble> >(numSEVars+numPBVars+2);
        for(int i = 0; i<numSEVars+numPBVars+2; i++)
        {
            k1[i] =  Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
            k2[i] =  Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
        }
  
        // // copy outarray to fp
        // for( int i = 0; i<nVariables; i++)
        // {
        //     Vmath::Vcopy(nSolutionPts, outay[i], 1, fp[i], 1);
        // }

        //Call RK-2 for PB and SE
        ComputeK(k1);

        
       //check if k1 is +ve:
       /* for( int i = 0; i<numSEVars+numPBVars+2; i++)
         {
            for(int j = 0; j<nScaledDomainPts; j++)
                if(k1[i][j]<0) {cout<<"\n k1 was negative after ComputeK for i = "<<i<<" and j = "<<j<<"val k1[i][j] = "<<k1[i][j];}
                }*/



        //go lower space -> higher space 
        
        // Scale all se vectors
        /*for( int i = 0; i<numSEVars; i++)
        {
            m_fields[0]->PhysInterp1DScaled(OneDptscaleS, se[i], seInterp[i]);
           
            for(int j = 0; j<nScaledDomainPts; j++)
                if(abs(seInterp[i][j])<1e-14)         
                    seInterp[i][j]=0; 
            
            // for(int j = 0; j< nScaledDomainPts; j++)
            // {   if(i !=5) 
            //         if(seInterp[i][j] < 0) 
            //         {
            //             cout<<"\n seinterp has -ve val at idx"<<
            //                 i<<" and "<<
            //                 j<<" val = "<<
            //                 seInterp[i][j];
            //             exit(0);  
            //      }
            //  }
            }

        // Scale all pb vectors
        for( int i = 0; i<numPBVars; i++)
        {
            m_fields[0]->PhysInterp1DScaled(OneDptscaleS, pb[i], pbInterp[i]);
             for(int j = 0; j< nScaledDomainPts; j++)
             {
                  if(abs(pbInterp[i][j])<1e-14)         
                    pbInterp[i][j]=0; 
            
                 // if(pbInterp[i][j]<0)
                 // {cout<<"\n pbinterp["<<i<<"]["<<j<<"] is -ve. Val = "<<pbInterp[i][j];exit(0);}
             }
        }*/
        

        // De-alias by sampling 1.5 times all fp vectors
        for( int i = 0; i<nVariables; i++)
        {
            m_fields[0]->PhysInterp1DScaled(OneDptscale, outarray[i], fpInterp[i]); 
            for(int j = 0; j<nSolutionPts; j++)
            {
                if(outarray[i][j]<0)
                {
                    cout<<"\n outarray["<<i<<"]["<<j<<"]="<<outarray[i][j];
                    //exit(0);
                }
            }
            for(int j = 0; j<nScaledDomainPts; j++)
            {
                 if(abs(fpInterp[i][j])<1e-14)
                     fpInterp[i][j] = 0; 
                 if(fpInterp[i][j]<0)
                 {
                     cout<<"\n fpInterp["<<i<<"]["<<j<<"] is -ve. Val = "<<fpInterp[i][j];
                     //exit(0); 
                 }
            }
            
        }        

        // Copy pb and se in tpbse
        Array<OneD, Array<OneD, NekDouble> > tpbse = Array<OneD, Array<OneD, NekDouble> >(numPBVars+numSEVars);
        // Copy PLPba and PLPsea in tpbapsea
        Array<OneD, Array<OneD, NekDouble> > tpbapsea = Array<OneD, Array<OneD, NekDouble> >(2);      
        int i;
        for(i = 0; i < numPBVars; i++)
        {
            tpbse[i] = Array<OneD, NekDouble>(nSolutionPts, 0.0);
            Vmath::Vcopy(nSolutionPts, pb[i], 1, tpbse[i], 1);
        }
        for(int j = 0; i < numSEVars + numPBVars; i++, j++)
        {
            tpbse[i] = Array<OneD, NekDouble>(nSolutionPts, 0.0);
           Vmath::Vcopy(nSolutionPts, se[j], 1, tpbse[i], 1);
        }
        tpbapsea[0] = Array<OneD, NekDouble>(nSolutionPts);
        tpbapsea[1] = Array<OneD, NekDouble>(nSolutionPts);
        
        Vmath::Vcopy(nSolutionPts, PLPba, 1, tpbapsea[0], 1);
        Vmath::Vcopy(nSolutionPts, PLPsea, 1, tpbapsea[1], 1);

            
        Array<OneD, NekDouble> temp = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
        Array<OneD, NekDouble> temp2 = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
        Array<OneD, NekDouble> Ae2 = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
        Array<OneD, NekDouble> Aadp = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);

       
        // Evaluation of reaction vectors in the higher space :
    
        // 0.Pmu = -kadh(Pmax-Psea)Pmu - (A1(e2) + A2([ADP]))Pmu 
        Vmath::Vmul(nScaledDomainPts, fpInterp[Pmu], 1, kadh, 1, fpReactInterp[Pmu], 1); // Pmu*kadh
        Vmath::Sadd(nScaledDomainPts, -1.0, PLPsea, 1, temp, 1); // Psea-1
        Vmath::Vmul(nScaledDomainPts, fpReactInterp[Pmu], 1, temp, 1, fpReactInterp[Pmu], 1); // -kadh(1-Psea)Pmu
        Vmath::Smul(nScaledDomainPts, Pmaxse, fpReactInterp[Pmu], 1, fpReactInterp[Pmu], 1); // -kadh(1-Psea)PmuPmax
  	
        //A1(e2) = [e2]/((1e-9/v2)+[e2]), v2=1.4e-6
        Vmath::Sadd(nScaledDomainPts, (1e-9/v2), fpInterp[FPe2], 1, Ae2, 1);
        Vmath::Vdiv(nScaledDomainPts, fpInterp[FPe2], 1, Ae2, 1, Ae2, 1);
        //A2([ADP]) = [adp]/((2e-6/va)+[ADP]), va = 5e-6
        Vmath::Sadd(nScaledDomainPts, (2e-6/va), fpInterp[ADP], 1, Aadp ,1); //kactadp[ADP]
        Vmath::Vdiv(nScaledDomainPts, fpInterp[ADP], 1, Aadp, 1, Aadp, 1); //((2e-6/va)+[e2])
        
        Vmath::Vmul(nScaledDomainPts, fpInterp[Pmu], 1, Ae2, 1, temp, 1);
        Vmath::Svtvp(nScaledDomainPts, -kact, temp, 1, fpReactInterp[Pmu], 1, fpReactInterp[Pmu], 1);
        Vmath::Vmul(nScaledDomainPts, fpInterp[Pmu], 1, Aadp, 1, temp, 1);
        Vmath::Svtvp(nScaledDomainPts, -kactadp, temp, 1, fpReactInterp[Pmu], 1, fpReactInterp[Pmu], 1);

        Vmath::Smul(nScaledDomainPts, tchar, fpReactInterp[Pmu], 1, fpReactInterp[Pmu], 1); //>
        /*cout<<"\nmax(Aadp) = "<<Vmath::Vmax(nScaledDomainPts, Aadp, 1)<<"\t max(Ae2) = "<<Vmath::Vmax(nScaledDomainPts, Ae2, 1)<<"\tmax(fpreactinterp[Pmu]="<<Vmath::Vmax(nScaledDomainPts, fpReactInterp[Pmu], 1);*/
        /*for(int i = 0; i<nScaledDomainPts ; i++) // fpReactInterp[Pmu] = max(0,fpReactInterp[Pmu])
        {
            if(fpReactInterp[Pmu][i]<0) 
            {
                cout<<"\n fpReactInterp[Pmu]["<<i<<"]="<<fpReactInterp[Pmu][i]; 
                fpReactInterp[Pmu][i] = 0; exit(0);
            }
            }*/
	
        // 1.Pma = -kadh(Pmax-Psea)Pma + A1(e2)Pmu +A2([ADP])Pmu - kcoh g(eta) PmaxPma   
        
        Vmath::Vmul(nScaledDomainPts, fpInterp[Pma]/*Pma*/, 1, kadh, 1, fpReactInterp[Pma], 1); // Pma*kadh
        Vmath::Sadd(nScaledDomainPts, -1.0, PLPsea, 1, temp, 1); // Psea-1 
        Vmath::Vmul(nScaledDomainPts, fpReactInterp[Pma], 1, temp, 1, fpReactInterp[Pma], 1); // -kadh(1-Psea)Pma 
        Vmath::Smul(nScaledDomainPts, Pmaxse, fpReactInterp[Pma], 1, fpReactInterp[Pma], 1); // -kadh(1-Psea)PmuPmax

        Vmath::Vmul(nScaledDomainPts, fpInterp[Pmu], 1, Ae2, 1, temp, 1);
        Vmath::Svtvp(nScaledDomainPts, kact, temp, 1, fpReactInterp[Pma], 1, fpReactInterp[Pma], 1);
        Vmath::Vmul(nScaledDomainPts, fpInterp[Pmu], 1, Aadp, 1, temp, 1);
        Vmath::Svtvp(nScaledDomainPts, kactadp, temp, 1, fpReactInterp[Pma], 1, fpReactInterp[Pma], 1);


        // g(eta) = max(0.0,betageta*(Peta(i,j)-etat)^3/(etastar^3 + (Peta(i,j)-etat)^3))
        NekDouble betageta = 1.0 / ( pow(1.0-etat,3) / (pow(etastar,3) + pow(1.0-etat,3)) );
	
        Vmath::Sadd(nScaledDomainPts, -etat, fpInterp[Peta], 1, temp, 1); //(Peta(i,j)-etat)
        Vmath::Vpow(nScaledDomainPts, temp, 1, 3.0, temp, 1); //(Peta(i,j)-etat)^3
        Vmath::Sadd(nScaledDomainPts, pow(etastar,3), temp, 1, temp2, 1);
        Vmath::Vdiv(nScaledDomainPts, temp, 1, temp2, 1, temp, 1); //temp = temp/temp2
        //temp = betageta*(Peta(i,j)-etat)^3/(etastar^3 + (Peta(i,j)-etat)^3))
        Vmath::Smul(nScaledDomainPts, betageta, temp, 1, geta, 1);
        for(int i = 0; i<nScaledDomainPts ; i++) // g(eta) = max(0,temp)
        {
            if(geta[i]<0) geta[i] = 0;
        }
        Vmath::Smul(nScaledDomainPts, -kcoh, fpInterp[Pma], 1, temp, 1); // temp2 = -Pma*kcoh
        Vmath::Vmul(nScaledDomainPts, geta, 1, temp, 1, temp, 1); // temp = -Pma*kcoh*g(eta)
        Vmath::Vadd(nScaledDomainPts, fpReactInterp[Pma], 1, temp, 1, fpReactInterp[Pma], 1); //done
        Vmath::Smul(nScaledDomainPts, tchar, fpReactInterp[Pma], 1, fpReactInterp[Pma], 1);//>
        /*for(int i = 0; i<nScaledDomainPts ; i++) // fpReactInterp[Pmu] = max(0,fpReactInterp[Pmu])
        {
            if(fpReactInterp[Pma][i]<0) 
            {
                cout<<"\n fpReactInterp[Pma]["<<i<<"]="<<fpReactInterp[Pma][i];
                fpReactInterp[Pma][i] = 0;exit(0);
            }
            }*/

        // 2.Peta = -gammaeta*Peta + gammaeta*((Pba+Psea))
        Vmath::Smul(nScaledDomainPts, -gammaeta, fpInterp[Peta], 1, fpReactInterp[Peta], 1);
        Vmath::Vadd(nScaledDomainPts, PLPba, 1, PLPsea, 1, temp, 1); //Pba+Psea
        Vmath::Smul(nScaledDomainPts, gammaeta, temp, 1, temp, 1); //gammaeta*((Pba+Psea))
        Vmath::Vadd(nScaledDomainPts, fpReactInterp[Peta], 1, temp, 1, fpReactInterp[Peta], 1); //done
        Vmath::Smul(nScaledDomainPts, tchar, fpReactInterp[Peta], 1, fpReactInterp[Peta], 1);//>
        /*for(int i = 0; i<nScaledDomainPts ; i++) // fpReactInterp[Pmu] = max(0,fpReactInterp[Pmu])
        {
            if(fpReactInterp[Peta][i]<0) 
            {
                cout<<"\n fpReactInterp[Peta]["<<i<<"]="<<fpReactInterp[Peta][i];exit(0);
                fpReactInterp[Peta][i] = 0;
            }
        }*/
        
        // 5. ADP = sigmaz/va - kadpinADP
        Vmath::Smul(nScaledDomainPts, 1/va, sigmaz, 1, fpReactInterp[ADP], 1);
        Vmath::Smul(nScaledDomainPts, kadpin, fpInterp[ADP], 1, temp, 1);
        Vmath::Vsub(nScaledDomainPts, fpReactInterp[ADP], 1, temp, 1, fpReactInterp[ADP], 1);
        Vmath::Smul(nScaledDomainPts, tchar, fpReactInterp[ADP], 1, fpReactInterp[ADP], 1);//>
        /*for(int i = 0; i<nScaledDomainPts ; i++) 
        {
            if(fpReactInterp[ADP][i]<0) 
            {
                
                cout<<"\n fpReactInterp[ADP]["<<i<<"]="<<fpReactInterp[ADP][i];exit(0);
                fpReactInterp[ADP][i] = 0;
            }
            }*/

        // 6. FPe2
        
        Vmath::Vadd(nScaledDomainPts, pb[PBe2], 1, pb[PBz5e2], 1, e2btotInterp, 1);
        Vmath::Svtvp(nScaledDomainPts, (s8/s2s), pb[PBz8e2], 1, e2btotInterp, 1, e2btotInterp, 1); // alpha*x+y
	
        Vmath::Smul(nScaledDomainPts, -ke2on*s2s, fpInterp[FPe2], 1, fpReactInterp[FPe2], 1);
        Vmath::Svtvm(nScaledDomainPts, N2starb*Pmaxb/s2s, PLPba, 1, e2btotInterp, 1, temp, 1); //alpha*x-y
        Vmath::Svtvp(nScaledDomainPts, N2starse*Pmaxse/s2s, PLPsea, 1, temp, 1, temp, 1); // alpha*x+y
        Vmath::Vmul(nScaledDomainPts, fpReactInterp[FPe2], 1, temp, 1, fpReactInterp[FPe2], 1);
        Vmath::Svtvp(nScaledDomainPts, ke2off*s2s/v2, pb[PBe2], 1, fpReactInterp[FPe2], 1, fpReactInterp[FPe2], 1); // alpha*x+y
        Vmath::Svtvp(nScaledDomainPts, -k2in, fpInterp[FPe2], 1, fpReactInterp[FPe2], 1, fpReactInterp[FPe2], 1); // alpha*x+y
        // (k14min+k14cat)*z8_e2*v7/v2
        Vmath::Svtvp(nScaledDomainPts, (kz8e2min + kz8e2cat)*(v8/v2), fpInterp[FPz8e2], 1, fpReactInterp[FPe2], 1, fpReactInterp[FPe2], 1);
        // (k12min+k12cat)*z5_e2*v5/v2
        Vmath::Svtvp(nScaledDomainPts, (kz5e2min + kz5e2cat)*(v5/v2), fpInterp[FPz5e2], 1, fpReactInterp[FPe2], 1, fpReactInterp[FPe2], 1);
        // (k18min+k18cat)*z7_e2*v7/v2
        Vmath::Svtvp(nScaledDomainPts, (kz7e2min + kz7e2cat)*(v7/v2), fpInterp[FPz7e2], 1, fpReactInterp[FPe2], 1, fpReactInterp[FPe2], 1);
        // - k18plu*z7*v7*e2
        Vmath::Smul(nScaledDomainPts, -(kz7e2plu)*v7, fpInterp[FPz7], 1, temp, 1); //-k18pluv7z7
        Vmath::Vvtvp(nScaledDomainPts, temp, 1, fpInterp[FPe2], 1, fpReactInterp[FPe2], 1, fpReactInterp[FPe2], 1);
        // - k12plu*z5*v5*e2
        Vmath::Smul(nScaledDomainPts, -(kz5e2plu)*v5, fpInterp[FPz5], 1, temp, 1); //-k12pluv5z5
        Vmath::Vvtvp(nScaledDomainPts, temp, 1, fpInterp[FPe2], 1, fpReactInterp[FPe2], 1, fpReactInterp[FPe2], 1);
        // - k14plu*z8*v8*e2
        Vmath::Smul(nScaledDomainPts, -(kz8e2plu)*v8, fpInterp[FPz8], 1, temp, 1); //-k14pluv8z8
        Vmath::Vvtvp(nScaledDomainPts, temp, 1, fpInterp[FPe2], 1, fpReactInterp[FPe2], 1, fpReactInterp[FPe2], 1);
        Vmath::Smul(nScaledDomainPts, tchar, fpReactInterp[FPe2], 1, fpReactInterp[FPe2], 1);//>
        /*for(int i = 0; i<nScaledDomainPts ; i++) 
         {
             if((fpReactInterp[FPe2][i])<0) 
             {
                 
                 cout<<"\n fpReactInterp[FPe2]["<<i<<"]="<<fpReactInterp[FPe2][i];exit(0);
                fpReactInterp[FPe2][i] = 0;
             }
             }*/
	
        // 7. FPz5
    
        Vmath::Vadd(nScaledDomainPts, pb[PBz5], 1, pb[PBz5e10], 1, ze5btotInterp, 1);
        Vmath::Svtvp(nScaledDomainPts, (s2/s5), pb[PBz5e2], 1, ze5btotInterp, 1, ze5btotInterp, 1); // alpha*x+y
        Vmath::Svtvp(nScaledDomainPts, (s10/s5), pb[PBz5e10], 1, ze5btotInterp, 1, ze5btotInterp, 1); // alpha*x+y
        Vmath::Vadd(nScaledDomainPts, pb[PBe5], 1, ze5btotInterp, 1, ze5btotInterp, 1);
        Vmath::Vadd(nScaledDomainPts, pb[PBapce5], 1, ze5btotInterp, 1, ze5btotInterp, 1);
        Vmath::Vadd(nScaledDomainPts, pb[PBpro], 1, ze5btotInterp, 1, ze5btotInterp, 1);
        Vmath::Vadd(nScaledDomainPts, pb[PBz2pro], 1, ze5btotInterp, 1, ze5btotInterp, 1);

        Vmath::Smul(nScaledDomainPts, -k5on*s5, fpInterp[FPz5], 1, fpReactInterp[FPz5], 1);
        Vmath::Svtsvtp(nScaledDomainPts, N5b*Pmaxb/s5, PLPba, 1, N5se*Pmaxse/s5, PLPsea, 1, temp, 1);

        Vmath::Vsub(nScaledDomainPts, temp, 1, ze5btotInterp, 1, temp, 1);

        Vmath::Vmul(nScaledDomainPts, fpReactInterp[FPz5], 1, temp, 1, fpReactInterp[FPz5], 1);

        Vmath::Svtvp(nScaledDomainPts, k5off*s5/v5, pb[PBz5], 1, fpReactInterp[FPz5], 1, fpReactInterp[FPz5], 1);

        Vmath::Smul(nScaledDomainPts, -(kz5e2plu)*v2, fpInterp[FPz5], 1, temp, 1); //- k12plu*z5*e2*v2 

        //cout<<" \n temp6 is : "<<Vmath::Vmax(nScaledDomainPts, temp, 1)<<" min = "<<Vmath::Vmin(nScaledDomainPts, temp, 1);
        Vmath::Vvtvp(nScaledDomainPts, temp, 1, fpInterp[FPe2], 1, fpReactInterp[FPz5], 1, fpReactInterp[FPz5], 1);

        Vmath::Svtvp(nScaledDomainPts, kz5e2min, fpInterp[FPz5e2], 1, fpReactInterp[FPz5], 1, fpReactInterp[FPz5], 1); //k12min*z5_e2

        // /+ N5*(Psea+Pba-oldp)/rk4tstep 
        Vmath::Vadd(nScaledDomainPts, PLPba, 1, PLPsea, 1, temp, 1);

        Vmath::Vsub(nScaledDomainPts, temp, 1, Pbold, 1, temp, 1);

        Vmath::Smul(nScaledDomainPts, N5/tchar, temp, 1, temp, 1);

        Vmath::Vadd(nScaledDomainPts, fpReactInterp[FPz5], 1, temp, 1, fpReactInterp[FPz5], 1);

        Vmath::Smul(nScaledDomainPts, tchar, fpReactInterp[FPz5], 1, fpReactInterp[FPz5], 1);//>

        /*for(int i = 0; i<nScaledDomainPts ; i++) 
        {
            if(fpReactInterp[FPz5][i]<0) 
            {
                
                cout<<"\n fpReactInterp[FPz5]["<<i<<"]="<<fpReactInterp[FPz5][i];exit(0);
                fpReactInterp[FPz5][i] = 0;
            }
            }*/

        // 8. FPz7
        Vmath::Svtvp(nScaledDomainPts, -kz7e10plu*v10, fpInterp[FPz7], 1, fpInterp[FPe10], 1, fpReactInterp[FPz7], 1); //-kz7e10plu v10 z7e10
        Vmath::Svtvp(nScaledDomainPts, kz7e10min, fpInterp[FPz7e10], 1, fpReactInterp[FPz7], 1, fpReactInterp[FPz7], 1); // +kz7e10min z7_e10
        Vmath::Svtvp(nScaledDomainPts, kz7e2min, fpInterp[FPz7e2], 1, fpReactInterp[FPz7], 1, fpReactInterp[FPz7], 1); // +kz7e2min z7_e2
        Vmath::Svtvp(nScaledDomainPts, -kz7e2plu*v2, fpInterp[FPz7], 1, fpInterp[FPe2], 1, temp, 1); // -kz7e2plu v2 z7e2
        Vmath::Vadd(nScaledDomainPts, temp, 1, fpReactInterp[FPz7], 1, fpReactInterp[FPz7], 1);
        Vmath::Smul(nScaledDomainPts, tchar, fpReactInterp[FPz7], 1, fpReactInterp[FPz7], 1);//>
        /*for(int i = 0; i<nScaledDomainPts ; i++) 
        {
            if(fpReactInterp[FPz7][i]<0) 
            {
                
                cout<<"\n fpReactInterp[FPz7]["<<i<<"]="<<fpReactInterp[FPz7][i];exit(0);
                fpReactInterp[FPz7][i] = 0;
            }
        }*/
	
        // 9. FPe10
        // (N10b Pba + N10se Psea - e10mtot - z10mtot)
        /*_________________ = ze10btot*/  	
        Vmath::Svtvp(nScaledDomainPts, (s8/s10), pb[PBz10tenstar], 1, pb[PBz10], 1, temp, 1); // z10ba + z10ba_tenstar*(s8/s10)
        Vmath::Svtvp(nScaledDomainPts, (s5/s10), pb[PBe10], 1, pb[PBpro], 1, temp2, 1); // e10ba + pro(s5/s10)
        Vmath::Vadd(nScaledDomainPts, temp, 1, temp2, 1, temp, 1);
        Vmath::Svtvp(nScaledDomainPts, (s8/s10), pb[PBz8e10], 1, temp, 1, temp, 1); // PBz8e10(s8/s10)
        Vmath::Svtvp(nScaledDomainPts, (s8/s10), pb[PBz10ten], 1, pb[PBz5e10], 1, temp2, 1); // pb[PBz5e10] + PBz10ten(s8/s10)
        Vmath::Vadd(nScaledDomainPts, temp, 1, temp2, 1, temp, 1);
        Vmath::Svtvp(nScaledDomainPts, (s5/s10), pb[PBz2pro], 1, temp, 1, temp, 1); // + pb[PBz2pro] (s5/s10)
	
        Vmath::Svtsvtp(nScaledDomainPts, N10b*Pmaxb/s10, PLPba, 1, N10se*Pmaxse/s10, PLPsea, 1, temp2, 1 );
        Vmath::Vsub(nScaledDomainPts, temp2, 1, temp, 1, temp, 1);// temp2 - temp
        Vmath::Svtvp(nScaledDomainPts, -k10on*s10, fpInterp[FPe10], 1, temp, 1, fpReactInterp[FPe10], 1);
        Vmath::Svtvp(nScaledDomainPts, -(s10/v10)*k10off, pb[PBe10], 1, fpReactInterp[FPe10], 1, fpReactInterp[FPe10], 1); // + k10offPBe10 s10/v10
        // -k10plu tfpi vf FPe10
        Vmath::Smul(nScaledDomainPts, -vf*ktfpie10plu, fpInterp[FPtfpi], 1, temp, 1);
        Vmath::Vvtvm (nScaledDomainPts, temp, 1, fpInterp[FPe10], 1, fpReactInterp[FPe10], 1, fpReactInterp[FPe10], 1); //vv+v
        // + k10min tfpi_e10 *(vf/v10)
        Vmath::Svtvp(nScaledDomainPts, (s10/v10)*ktfpie10min, fpInterp[FPtfpie10], 1, fpReactInterp[FPe10], 1, fpReactInterp[FPe10], 1);
        // + (k1min+k1cat)*z7_e10 *v7/v10
        Vmath::Svtvp(nScaledDomainPts, (v7/v10)*(kz7e10min + kz7e10cat), fpInterp[FPz7e10], 1, fpReactInterp[FPe10], 1, fpReactInterp[FPe10], 1);
        // -k1plu z7 v7 e10
        Vmath::Smul(nScaledDomainPts, -kz7e10plu*v7, fpInterp[FPz7], 1, temp, 1);
        Vmath::Vvtvm (nScaledDomainPts, temp, 1, fpInterp[FPe10], 1, fpReactInterp[FPe10], 1, fpReactInterp[FPe10], 1);
        // -k10in e10
        Vmath::Svtvp(nScaledDomainPts, -k10in, fpInterp[FPe10], 1, fpReactInterp[FPe10], 1, fpReactInterp[FPe10], 1);
        Vmath::Smul(nScaledDomainPts, tchar, fpReactInterp[FPe10], 1, fpReactInterp[FPe10], 1);//>
	    /*for(int i = 0; i<nScaledDomainPts ; i++) 
           {
               if(fpReactInterp[FPe10][i]<0)
               {
                
                   cout<<"\n fpReactInterp[FPe10]["<<i<<"]="<<fpReactInterp[FPe10][i];exit(0);
                   fpReactInterp[FPe10][i] = 0;
               }
               }*/
       
        // 10. FPz8
        //-k8on z8 s8 (N8b Pmaxb Pba (1/s8) + N8sePmaxse Psea (1/s8) -ze8btot)
        //ze8btot = 
        //PBz8 + PBz8e10 + pb[PBz8e2] + PBe8 + PBapce8 + PBten + PBz10ten + PBtenstar + pb[PBz10tenstar]
        Vmath::Vadd(nScaledDomainPts, pb[PBz8], 1, pb[PBz8e10], 1, ze8btotInterp, 1);
        Vmath::Vadd(nScaledDomainPts, ze8btotInterp, 1, pb[PBz8e2], 1, ze8btotInterp, 1);
        Vmath::Vadd(nScaledDomainPts, ze8btotInterp, 1, pb[PBe8], 1, ze8btotInterp, 1);
        Vmath::Vadd(nScaledDomainPts, ze8btotInterp, 1, pb[PBapce8], 1, ze8btotInterp, 1);
        Vmath::Vadd(nScaledDomainPts, ze8btotInterp, 1, pb[PBten], 1, ze8btotInterp, 1);
        Vmath::Vadd(nScaledDomainPts, ze8btotInterp, 1, pb[PBz10ten], 1, ze8btotInterp, 1);
        Vmath::Vadd(nScaledDomainPts, ze8btotInterp, 1, pb[PBtenstar], 1, ze8btotInterp, 1);
        Vmath::Vadd(nScaledDomainPts, ze8btotInterp, 1, pb[PBz10tenstar], 1, ze8btotInterp, 1);
    
        Vmath::Svtsvtp(nScaledDomainPts, N8b*Pmaxb/s8, PLPba, 1, N8se*Pmaxse/s8, PLPsea, 1, fpReactInterp[FPz8], 1);
        Vmath::Vsub(nScaledDomainPts, fpReactInterp[FPz8], 1, ze8btotInterp, 1, fpReactInterp[FPz8], 1);
        Vmath::Smul(nScaledDomainPts, -k8on*s8, fpInterp[FPz8], 1, temp, 1);
        Vmath::Vmul(nScaledDomainPts, temp, 1, fpReactInterp[FPz8], 1, fpReactInterp[FPz8], 1);
        // + k8off PBz8 (s8/v8)
        Vmath::Svtvp(nScaledDomainPts, k8off*(s8/v8), pb[PBz8], 1, fpReactInterp[FPz8], 1, fpReactInterp[FPz8], 1);
        // -kz8e2plu z8 e2 v2
        Vmath::Smul(nScaledDomainPts, -kz8e2plu*v2, fpInterp[FPz8], 1, temp, 1);
        Vmath::Vmul(nScaledDomainPts, temp, 1, fpInterp[FPe2], 1, temp, 1);
        Vmath::Vadd(nScaledDomainPts, fpReactInterp[FPz8], 1, temp, 1, fpReactInterp[FPz8], 1);
        // + kz8e2min z8_e2
        Vmath::Svtvp(nScaledDomainPts, kz8e2min, fpInterp[FPz8e2], 1, fpReactInterp[FPz8], 1, fpReactInterp[FPz8],1);
        Vmath::Smul(nScaledDomainPts, tchar, fpReactInterp[FPz8], 1, fpReactInterp[FPz8], 1);//>
        /*for(int i = 0; i<nScaledDomainPts ; i++) 
        {
            if(fpReactInterp[FPz8][i]<0) 
            {
                
                cout<<"\n fpReactInterp[FPz8]["<<i<<"]="<<fpReactInterp[FPz8][i];exit(0);
                fpReactInterp[FPz8][i] = 0;
            }
            }*/

        // 11. FPz5e2
        // kz5e2plu*z5*e2*v2 
        Vmath::Smul(nScaledDomainPts, kz5e2plu*v2, fpInterp[FPz5], 1, temp, 1);
        Vmath::Vmul(nScaledDomainPts, temp, 1, fpInterp[FPe2], 1, fpReactInterp[FPz5e2], 1);
        // - (kz5e2min + kz5e2cat)z5_e2
        Vmath::Svtvp(nScaledDomainPts, -(kz5e2min+kz5e2cat), fpInterp[FPz5e2], 1, fpReactInterp[FPz5e2], 1, fpReactInterp[FPz5e2],1);
        Vmath::Smul(nScaledDomainPts, tchar, fpReactInterp[FPz5e2], 1, fpReactInterp[FPz5e2], 1);//>
         for(int i = 0; i<nScaledDomainPts ; i++) 
        {
            if(fpReactInterp[FPz5e2][i]<0) fpReactInterp[FPz5e2][i] = 0;
        }
    
        // 12. FPz7e2
        // kz7e2plu z7 e2 v2
        Vmath::Smul(nScaledDomainPts, kz7e2plu*v2, fpInterp[FPz7], 1, temp, 1);
        Vmath::Vmul(nScaledDomainPts, temp, 1, fpInterp[FPe2], 1, fpReactInterp[FPz7e2], 1);
        // -(kz7e2min + kz7e2cat)z7_e2
        Vmath::Svtvp(nScaledDomainPts, -(kz7e2min+kz7e2cat), fpInterp[FPz7e2], 1, fpReactInterp[FPz7e2], 1, fpReactInterp[FPz7e2],1);
        Vmath::Smul(nScaledDomainPts, tchar, fpReactInterp[FPz7e2], 1, fpReactInterp[FPz7e2], 1);//>
        /*for(int i = 0; i<nScaledDomainPts ; i++) 
        {
            if(fpReactInterp[FPz7e2][i]<0) fpReactInterp[FPz7e2][i] = 0;
            }*/

        //13. FPz7e10
        // kz7e10plu z7 e10 v10
        Vmath::Smul(nScaledDomainPts, kz7e10plu*v10, fpInterp[FPz7], 1, temp, 1);
        Vmath::Vmul(nScaledDomainPts, temp, 1, fpInterp[FPe10], 1, fpReactInterp[FPz7e10], 1);
        // -(kz7e10min + kz7e10cat)z7_e10
        Vmath::Svtvp(nScaledDomainPts, -(kz7e10min+kz7e10cat), fpInterp[FPz7e10], 1, fpReactInterp[FPz7e10], 1, fpReactInterp[FPz7e10],1);
        Vmath::Smul(nScaledDomainPts, tchar, fpReactInterp[FPz7e10], 1, fpReactInterp[FPz7e10], 1);//>
        /*for(int i = 0; i<nScaledDomainPts ; i++) 
        {
            if(fpReactInterp[FPz7e10][i]<0) fpReactInterp[FPz7e10][i] = 0;
            }*/

        //14. FPz8e2
        // kz8e2plu z8 e2 v2
        Vmath::Smul(nScaledDomainPts, kz8e2plu*v2, fpInterp[FPz8], 1, temp, 1);
        Vmath::Vmul(nScaledDomainPts, temp, 1, fpInterp[FPe2], 1, fpReactInterp[FPz8e2], 1);
        // -(kz8e2min + kz8e2cat)z8_e2
        Vmath::Svtvp(nScaledDomainPts, -(kz8e2min+kz8e2cat), fpInterp[FPz8e2], 1, fpReactInterp[FPz8e2], 1, fpReactInterp[FPz8e2],1);
        Vmath::Smul(nScaledDomainPts, tchar, fpReactInterp[FPz8e2], 1, fpReactInterp[FPz8e2], 1);//>
        /*for(int i = 0; i<nScaledDomainPts ; i++) 
        {
            if(fpReactInterp[FPz8e2][i]<0) fpReactInterp[FPz8e2][i] = 0;
            }*/

        //15. FPtfpi
        // -ktfpie10plu tfpi e10 v10
        Vmath::Smul(nScaledDomainPts, -ktfpie10plu*v2, fpInterp[FPtfpi], 1, temp, 1);
        Vmath::Vmul(nScaledDomainPts, temp, 1, fpInterp[FPe10], 1, fpReactInterp[FPtfpi], 1);
        // + ktfpie10min tfpi_e10
        Vmath::Svtvp(nScaledDomainPts, ktfpie10min, fpInterp[FPtfpie10], 1, fpReactInterp[FPtfpi], 1, fpReactInterp[FPtfpi],1);
        Vmath::Smul(nScaledDomainPts, tchar, fpReactInterp[FPtfpi], 1, fpReactInterp[FPtfpi], 1);//>
        for(int i = 0; i<nScaledDomainPts ; i++) 
        {
            if(fpReactInterp[FPtfpi][i]<0) fpReactInterp[FPtfpi][i] = 0;
        }

        //16. FPtfpie10
        //ktfpie10plu tfpi e10 v10
        Vmath::Smul(nScaledDomainPts, ktfpie10plu*v10, fpInterp[FPtfpi], 1, temp, 1);
        Vmath::Vmul(nScaledDomainPts, temp, 1, fpInterp[FPe10], 1, fpReactInterp[FPtfpie10], 1);    
        //-ktfpie10min tfpi_e10
        Vmath::Svtvp(nScaledDomainPts, -ktfpie10min, fpInterp[FPtfpie10], 1, fpReactInterp[FPtfpie10], 1, fpReactInterp[FPtfpie10],1);
        Vmath::Smul(nScaledDomainPts, tchar, fpReactInterp[FPtfpie10], 1, fpReactInterp[FPtfpie10], 1);//>
        /*for(int i = 0; i<nScaledDomainPts ; i++) 
        {
            if(fpReactInterp[FPtfpie10][i]<0) fpReactInterp[FPtfpie10][i] = 0;
            }*/
        
        //17. FPz2
        Vmath::Svtvp(nScaledDomainPts, (s5/s2), pb[PBz2pro], 1, pb[PBz2], 1, z2btotInterp, 1);
        //N2b*Pmaxb*Pba/s2 + N2se*Pmaxse*Psea/s2-z2btot
        Vmath::Svtsvtp(nScaledDomainPts, N2b*Pmaxb/s2, PLPba, 1, N2se*Pmaxse/s2, PLPsea, 1, temp, 1);
        Vmath::Vsub(nScaledDomainPts, temp, 1, z2btotInterp, 1, temp, 1);
        //-k2on*z2*s2*temp
        Vmath::Smul(nScaledDomainPts, -k2on*s2, fpInterp[FPz2], 1, temp2, 1);
        Vmath::Vmul(nScaledDomainPts, temp, 1, temp2, 1, temp, 1);
        //+k2off*PBz2*s2/v2
        Vmath::Svtvp(nScaledDomainPts, k2off*s2/v2, pb[PBz2], 1, temp, 1, temp, 1);
        Vmath::Smul(nScaledDomainPts, tchar, temp, 1, fpReactInterp[FPz2], 1);//>
        /*for(int i = 0; i<nScaledDomainPts ; i++) 
        {
            if(fpReactInterp[FPz2][i]<0) fpReactInterp[FPz2][i] = 0;
            }*/

        //18. FPapc
        //-kapce5pbplu*FPapc*pb[PBe5]*s5 
        Vmath::Vmul(nScaledDomainPts, fpInterp[FPapc], 1, pb[PBe5], 1, fpReactInterp[FPapc], 1);
        Vmath::Smul(nScaledDomainPts, -kapce5pbplu*s5, fpReactInterp[FPapc], 1, fpReactInterp[FPapc], 1);
        //+(kapce5pbmin+kapce5pbcat)*PBapce5*s5/v2
        Vmath::Svtvp(nScaledDomainPts, (kapce5pbmin+kapce5pbcat)*s5/v2, pb[PBapce5], 1, fpReactInterp[FPapc], 1, fpReactInterp[FPapc], 1);
        //-kapce8pbplu*FPapc*PBe8*s8 
        Vmath::Vmul(nScaledDomainPts, fpInterp[FPapc], 1, pb[PBe8], 1, temp2, 1);
        Vmath::Svtvp(nScaledDomainPts, -kapce8pbplu*s8, temp2, 1, fpReactInterp[FPapc], 1, fpReactInterp[FPapc], 1);
        //+(kapce8pbmin+kapce8pbcat)*PBapce8*s8/v2
        Vmath::Svtvp(nScaledDomainPts, (kapce8pbmin+kapce8pbcat)*s8/v2, pb[PBapce8], 1, fpReactInterp[FPapc], 1, fpReactInterp[FPapc], 1);
        Vmath::Smul(nScaledDomainPts, tchar, temp, 1, fpReactInterp[FPapc], 1);//>
        /*for(int i = 0; i<nScaledDomainPts ; i++) 
        {
            if(fpReactInterp[FPapc][i]<0) fpReactInterp[FPapc][i] = 0;
            }*/

        //19. FPz10
        //N10b*Pmaxb*Pba/s10 + N10se*Pmaxse*Psea/s10
        Vmath::Svtsvtp(nScaledDomainPts, N10b*Pmaxb/s10, PLPba, 1, N10se*Pmaxse/s10, PLPsea, 1, fpReactInterp[FPz10], 1);
        //-ze10btot
        Vmath::Vsub(nScaledDomainPts, fpReactInterp[FPz10], 1, ze10btotInterp, 1, fpReactInterp[FPz10], 1);
        //-k10on*z10*s10*temp
        Vmath::Vmul(nScaledDomainPts, fpInterp[FPz10], 1, fpReactInterp[FPz10], 1, fpReactInterp[FPz10], 1);
        Vmath::Smul(nScaledDomainPts, -k10on*s10, fpReactInterp[FPz10], 1, fpReactInterp[FPz10], 1);
        //+k10zoff*PBz10*s10/v10
        Vmath::Svtvp(nScaledDomainPts, (k10zoff)*s10/v10, pb[PBz10], 1, fpReactInterp[FPz10], 1, fpReactInterp[FPz10], 1);
        Vmath::Smul(nScaledDomainPts, tchar, fpReactInterp[FPz10], 1, fpReactInterp[FPz10], 1);//>
        /*for(int i = 0; i<nScaledDomainPts ; i++) 
        {
            if(fpReactInterp[FPz10][i]<0) fpReactInterp[FPz10][i] = 0;
            }*/

        //20. FPe5
        //N5b*Pmaxb*Pba/s5 + N5se*Pmaxse*Psea/s5
        Vmath::Svtsvtp(nScaledDomainPts, N5b*Pmaxb/s5, PLPba, 1, N5se*Pmaxse/s5, PLPsea, 1, fpReactInterp[FPe5], 1);
        //-ze5btot
        Vmath::Vsub(nScaledDomainPts, fpReactInterp[FPe5], 1, ze5btotInterp, 1, fpReactInterp[FPe5], 1);
        //-k5on*FPe5*s5*temp
        Vmath::Vmul(nScaledDomainPts, fpInterp[FPe5], 1, fpReactInterp[FPe5], 1, fpReactInterp[FPe5], 1);
        Vmath::Smul(nScaledDomainPts, -k5on*s5, fpReactInterp[FPe5], 1, fpReactInterp[FPe5], 1);
        //+k5off*pb[PBe5]*s5/v5
        Vmath::Svtvp(nScaledDomainPts, (k5off)*s5/v5, pb[PBe5], 1, fpReactInterp[FPe5], 1, fpReactInterp[FPe5], 1);
        //+kz5e2cat*FPz5e2
        Vmath::Svtvp(nScaledDomainPts, kz5e2cat, fpInterp[FPz5e2], 1, fpReactInterp[FPe5], 1, fpReactInterp[FPe5], 1);
        Vmath::Smul(nScaledDomainPts, tchar, fpReactInterp[FPe5], 1, fpReactInterp[FPe5], 1);//>
        /*for(int i = 0; i<nScaledDomainPts ; i++) 
        {
            if(fpReactInterp[FPe5][i]<0) fpReactInterp[FPe5][i] = 0;
            }*/
        
        //22. FPe8
        //N8b*Pmaxb*Pba/s8 + N8se*Pmaxse*Psea/s8
        Vmath::Svtsvtp(nScaledDomainPts, N8b*Pmaxb/s8, PLPba, 1, N8se*Pmaxse/s8, PLPsea, 1, fpReactInterp[FPe8], 1);
        //-ze8btotInterp
        Vmath::Vsub(nScaledDomainPts, fpReactInterp[FPe8], 1, ze8btotInterp, 1, fpReactInterp[FPe8], 1);
        //-k8on*FPe8*s8*temp
        Vmath::Vmul(nScaledDomainPts, fpInterp[FPe8], 1, fpReactInterp[FPe8], 1, fpReactInterp[FPe8], 1);
        Vmath::Smul(nScaledDomainPts, -k8on*s8, fpReactInterp[FPe8], 1, fpReactInterp[FPe8], 1);
        //+k8off*PBe8*s8/v8 
        Vmath::Svtvp(nScaledDomainPts, k8off*s8/v8, pb[PBe8], 1, fpReactInterp[FPe8], 1, fpReactInterp[FPe8], 1);
        //+k14cat*FPz8e2
        Vmath::Svtvp(nScaledDomainPts, kz8e2cat, fpInterp[FPz8e2], 1, fpReactInterp[FPe8], 1, fpReactInterp[FPe8], 1);
        Vmath::Smul(nScaledDomainPts, tchar, fpReactInterp[FPe8], 1, fpReactInterp[FPe8], 1);//>
        /*for(int i = 0; i<nScaledDomainPts ; i++) 
        {
            if(fpReactInterp[FPe8][i]<0) fpReactInterp[FPe8][i] = 0;
            }*/

        //23. FPe9
        //ze9btot: PBz9 + PBe9
        Vmath::Vmul(nScaledDomainPts, pb[PBe9], 1, pb[PBz9], 1, ze9btotInterp, 1);
        //+PBz10ten*(s8/s9) 
        Vmath::Svtvp(nScaledDomainPts, (s8/s9), pb[PBz10ten], 1, ze9btotInterp, 1, ze9btotInterp, 1);
        //+PBten*(s8/s9)
        Vmath::Svtvp(nScaledDomainPts, (s8/s9), pb[PBten], 1, ze9btotInterp, 1, ze9btotInterp, 1);
        
        //N9b*Pmaxb*Pba/s9 + N9se*Pmaxse*Psea/s9
        Vmath::Svtsvtp(nScaledDomainPts, N9b*Pmaxb/s9, PLPba, 1, N9se*Pmaxse/s9, PLPsea, 1, fpReactInterp[FPe9], 1);
        //-ze9btotInterp
        Vmath::Vsub(nScaledDomainPts, fpReactInterp[FPe9], 1, ze9btotInterp, 1, fpReactInterp[FPe9], 1);
        //-k9on*FPe9*s9*temp
        Vmath::Vmul(nScaledDomainPts, fpInterp[FPe9], 1, fpReactInterp[FPe9], 1, fpReactInterp[FPe9], 1);
        Vmath::Smul(nScaledDomainPts, -k9on*s9, fpReactInterp[FPe9], 1, fpReactInterp[FPe9], 1);
        //+k9off*PBe9*s9/v9 
        Vmath::Svtvp(nScaledDomainPts, k9off*s9/v9, pb[PBe9], 1, fpReactInterp[FPe9], 1, fpReactInterp[FPe9], 1);
        //-k9in*FPe9
        Vmath::Svtvp(nScaledDomainPts, -k9in, fpInterp[FPe9], 1, fpReactInterp[FPe9], 1, fpReactInterp[FPe9], 1);       
        //N9starb*Pmaxb*Pba + N9starse*Pmaxse*Psea            
        Vmath::Svtsvtp(nScaledDomainPts, N9starb*Pmaxb, PLPba, 1, N9starse*Pmaxse, PLPsea, 1, temp, 1);
        //- PBe9star*s91 
        Vmath::Svtvp(nScaledDomainPts, -s91, pb[PBe9star], 1, temp, 1, temp, 1);
        //-PBtenstar*s8 
        Vmath::Svtvp(nScaledDomainPts, -s8, pb[PBtenstar], 1, temp, 1, temp, 1); 
        //- pb[PBz10tenstar]*s8
        Vmath::Svtvp(nScaledDomainPts, -s8, pb[PBz10tenstar], 1, temp, 1, temp, 1);   
        //-k9on*FPe9*temp
        Vmath::Vmul(nScaledDomainPts, temp, 1, fpInterp[FPe9], 1, temp, 1);
        Vmath::Svtvp(nScaledDomainPts, -k9on, temp, 1, fpReactInterp[FPe9], 1, fpReactInterp[FPe9], 1);
        //+k9off*PBe9star*s91/v9
        Vmath::Svtvp(nScaledDomainPts, k9off*s91/v9, pb[PBe9star], 1, fpReactInterp[FPe9], 1, fpReactInterp[FPe9], 1);
        Vmath::Smul(nScaledDomainPts, tchar, fpReactInterp[FPe9], 1, fpReactInterp[FPe9], 1);//>                                                                                                     
        /*for(int i = 0; i<nScaledDomainPts ; i++) 
        {
            if(fpReactInterp[FPe9][i]<0) fpReactInterp[FPe9][i] = 0;
            }*/

        //24. FPz9
        //temp = N9b*Pmaxb*Pba/s9 + N9se*Pmaxse*Psea/s9
        Vmath::Svtsvtp(nScaledDomainPts, N9b*Pmaxb/s9, PLPba, 1, N9se*Pmaxse/s9, PLPsea, 1, fpReactInterp[FPz9], 1);
        //-ze9btotInterp
        Vmath::Vsub(nScaledDomainPts, fpReactInterp[FPz9], 1, ze9btotInterp, 1, fpReactInterp[FPz9], 1);
        //-k9on*z9*s9*temp
        Vmath::Vmul(nScaledDomainPts, fpInterp[FPz9], 1, fpReactInterp[FPz9], 1, fpReactInterp[FPz9], 1);
        Vmath::Smul(nScaledDomainPts, -k9on*s9, fpReactInterp[FPz9], 1, fpReactInterp[FPz9], 1);
        //+k9off*PBz9*s9/v9
        Vmath::Svtvp(nScaledDomainPts, k9off*s9/v9, pb[PBz9], 1, fpReactInterp[FPz9], 1, fpReactInterp[FPz9], 1);
        Vmath::Smul(nScaledDomainPts, tchar, fpReactInterp[FPz9], 1, fpReactInterp[FPz9], 1);//>  
        /*for(int i = 0; i<nScaledDomainPts ; i++) 
        {
            if(fpReactInterp[FPz9][i]<0) fpReactInterp[FPz9][i] = 0;
            }*/

        //25. FPe7
        //k1cat*z7_e10 + k18cat*z7_e2
        Vmath::Svtsvtp(nScaledDomainPts, kz7e10cat, fpInterp[FPz7e10], 1, kz7e2cat, fpInterp[FPz7e2], 1, fpReactInterp[FPe7], 1);
        Vmath::Smul(nScaledDomainPts, tchar, fpReactInterp[FPe7], 1, fpReactInterp[FPe7], 1);//>  
        /*for(int i = 0; i<nScaledDomainPts ; i++) 
        {
            if(fpReactInterp[FPz9][i]<0) fpReactInterp[FPz9][i] = 0;
            }*/
       
        // update Cvals[$] with Cvals[$] + dt k1[$]
        for(i = 0; i < numPBVars; i++)
        {
            Vmath::Svtvp(nSolutionPts, tchar, k1[i], 1, 
                         pb[i], 1, pb[i], 1);
        }
    
        for(int j = 0; i < numSEVars + numPBVars; i++, j++)
        {
            Vmath::Svtvp(nSolutionPts, tchar, k1[i], 1, 
                         se[j], 1, se[j], 1);
        }
    
        Vmath::Svtvp(nSolutionPts, tchar, k1[i++], 1, 
                     PLPba, 1, PLPba, 1);
        Vmath::Svtvp(nSolutionPts, tchar, k1[i++], 1, 
                     PLPsea, 1, PLPsea, 1);

        //cout<<"\n Before computek(k2):";
        //check if pbinterp is all +ve:

        // update fp to fp+reactInterp at current dt
        // 1. fpInterp = fpInterp+fpReactInterp
        for (int i = 0; i < nVariables; ++i) 
        {
            m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, 
                                                        fpReactInterp[i], 
                                                        outArrayReact[i]);       
        }

        //add react term to outarray
        for (int i = 0; i < nVariables; ++i)
        {
            Vmath::Vadd(nSolutionPts, &outarray[i][0], 1, 
                        &outArrayReact[i][0], 1, &outarray[i][0], 1);            
             for(int k = 0; k<nSolutionPts; k++)
             {
                 if(abs(outarray[i][k])<1e-14) 
                 {
                     outarray[i][k] = 0;
                 }
              
                 if(outarray[i][k]<0) 
                 {
                     //outarray[i][k] = 0;
                     cout<<"\noutArray negative after adding react term! for i = "<<i<<
                         " and k = "<<k<<" val = "<<
                         outarray[i][k];
                     //exit(0);
                 }            
             }
        }

        
        // // 2. fp = outarray
        // for (int i = 0; i < nVariables; ++i) 
        // {
        //     Vmath::Vcopy(nSolutionPts, outarray[i], 1, fp[i], 1);
        // }

        ComputeK(k2);

        //SSPRK-2 : u^{2} = 0.5 u^{0} + 0.5(u^{1}) + 0.5dtf(u^{1})
        //currently, pb and se has u^{1} & tpbse, etc 
        //has u^{0} & k2 has f(u^1)
        
        for(i = 0; i < numPBVars; i++)
        {
            Vmath::Svtsvtp(nSolutionPts, 0.5, pb[i], 1, 
                           0.5, tpbse[i], 1, pb[i], 1);
            Vmath::Svtvp(nSolutionPts, 0.5*tchar, k2[i], 1, 
                         pb[i], 1, pb[i], 1); 
        }

        for(int j = 0; i < numSEVars + numPBVars; i++, j++)
        {
            Vmath::Svtsvtp(nSolutionPts, 0.5, se[j], 1, 
                           0.5, tpbse[i], 1, se[j], 1);
 
            Vmath::Svtvp(nSolutionPts, 0.5*tchar, k2[i], 1, 
                         se[j], 1, se[j], 1); 
                           
        }
        
        Vmath::Svtsvtp(nSolutionPts, 0.5, PLPba, 1, 
                       0.5, tpbapsea[0], 1, PLPba, 1);
        Vmath::Svtvp(nSolutionPts, 0.5*tchar, k2[i++], 1, 
                     PLPba, 1, PLPba, 1); 
      
        Vmath::Svtsvtp(nSolutionPts, 0.5, PLPsea, 1, 
                       0.5, tpbapsea[1], 1, PLPsea, 1);
        Vmath::Svtvp(nSolutionPts, 0.5*tchar, k2[i], 1, 
                     PLPsea, 1, PLPsea, 1); 

        // //RK-2 Heun's
        // // Copy back tpbse and tpbasea to pb, se and PLPba and PLPsea
        // for( i = 0; i< numPBVars; i++)
        // {
        //     Vmath::Vcopy(nSolutionPts, tpbse[i], 1, pb2[i], 1);
        // }
        // for ( int j = 0; i<numPBVars+numSEVars; i++, j++)
        // {
        //     Vmath::Vcopy(nSolutionPts, tpbse[i], 1, se2[j], 1);
        // }
            
        // // k1 = dt * k1 and k2 = dt * k2
        // for(i = 0; i < numPBVars+numSEVars+2; i++)
        // {
        //     Vmath::Smul(nSolutionPts, tchar, k1[i], 1, k1[i], 1);
        //     Vmath::Smul(nSolutionPts, tchar, k2[i], 1, k2[i], 1);
        // }        
        // //apply RK2 formula that uses k1 and k2:
        // //(a) k1 = k1 + k2
        // for(i = 0; i < numPBVars+numSEVars+2; i++)
        // {
        //     Vmath::Vadd(nSolutionPts, k1[i], 1, k2[i], 1, k1[i], 1);
        // }
        // // (b) (1/2)*(k1+k2), store in k1
        // for(i = 0; i < numPBVars+numSEVars+2; i++)
        // {
        //     Vmath::Smul(nSolutionPts, 0.5, k1[i], 1, k1[i], 1);
        // }
        // // (c) update Cvals[$] with  Cvals[$] + k1
        // for(i = 0; i < numPBVars; i++)
        // {
        //     Vmath::Vadd(nSolutionPts, pb2[i], 1, k1[i], 1, pb2[i], 1);
        // }
        // for(int j = 0; i < numSEVars + numPBVars; i++, j++)
        // {
        //     Vmath::Vadd(nSolutionPts, k1[i], 1, se2[j], 1, se2[j], 1);
        // }

        //RK-2 Ralston's
        // Copy back tpbse and tpbasea to pb, se and PLPba and PLPsea
        for( i = 0; i< numPBVars; i++)
        {
            Vmath::Vcopy(nSolutionPts, tpbse[i], 1, pb2[i], 1);
        }
        for ( int j = 0; i<numPBVars+numSEVars; i++, j++)
        {
            Vmath::Vcopy(nSolutionPts, tpbse[i], 1, se2[j], 1);
        }
            
        // k1 = dt * k1 and k2 = dt * k2
        for(i = 0; i < numPBVars+numSEVars+2; i++)
        {
            Vmath::Smul(nSolutionPts, tchar, k1[i], 1, k1[i], 1);
            Vmath::Smul(nSolutionPts, tchar, k2[i], 1, k2[i], 1);
        }        
        //apply RK2 formula that uses k1 and k2:
        //(a) k1 = (1/4)k1 + (3/4)k2
        for(i = 0; i < numPBVars+numSEVars+2; i++)
        {
            //Vmath::Svtsvtp(nSolutionPts, 1/4, k1[i], 1, 3/4, k2[i], 1, k1[i], 1);
            Vmath::Smul(nSolutionPts, 0.25, k1[i], 1, k1[i], 1);
            Vmath::Smul(nSolutionPts, 0.75, k2[i], 1, k2[i], 1);
            Vmath::Vadd(nSolutionPts, k1[i], 1, k2[i], 1, k1[i], 1);
        }
        // (c) update Cvals[$] with  Cvals[$] + k1
        for(i = 0; i < numPBVars; i++)
        {
            Vmath::Vadd(nSolutionPts, pb2[i], 1, k1[i], 1, pb2[i], 1);
        }
        for(int j = 0; i < numSEVars + numPBVars; i++, j++)
        {
            Vmath::Vadd(nSolutionPts, k1[i], 1, se2[j], 1, se2[j], 1);
        }


        // check if ssp-rk2 is same as rk-2
        for (int i = 0; i < numSEVars; ++i) 
        {
            for(int j = 0; j<nSolutionPts; j++)
            {
                if( se2[i][j]<0)//se2[i][j]-se[i][j]>1e-14 ||
                {
                    cout<<"se["<<i<<"]["<<j<<"]="<<se[i][j]<<" and se2["<<i<<"]["<<j<<"] = "<<se2[i][j]<< " diff = "<<se[i][j]-se2[i][j];
                    exit(0);
                }
            }
        }
        for (int i = 0; i < numPBVars; ++i) 
        {
            for(int j = 0; j<nSolutionPts; j++)
            {
                if( pb2[i][j]<0)//pb2[i][j]-pb[i][j]>1e-14||
                {
                    cout<<"pb["<<i<<"]["<<j<<"]="<<pb[i][j]<<" and pb2["<<i<<"]["<<j<<"] = "<<pb2[i][j];
                    exit(0);
                }
            }
        }

        for (int i = 0; i < numSEVars; ++i) 
        {
            for(int j = 0; j<nSolutionPts; j++)
            {
                if(abs(se[i][j])<1e-14)
                    se[i][j] = 0;
                if(se[i][j]<0) {cout<<"\n1se["<<i<<"]["<<j<<"]="
                                          <<se[i][j];  exit(0); }
                 
            }
        }
        for (int i = 0; i < numPBVars; ++i) 
        {
            
            for(int j = 0; j<nSolutionPts; j++)
            {

                if(abs(pb[i][j])<1e-14)
                    pb[i][j] = 0;
                if(pb[i][j]<0) 
                {
                    cout<<"\n1pb["<<i<<"]["<<j<<"]="<<pb[i][j];  
                    exit(0);
                }
        
            }
        }
 
        //Update PboldInterp here for next timestep
        //m_fields[0]->PhysInterp1DScaled(OneDptscale, Pbold, PboldInterp);

        if (prev == m_session->GetParameter("IO_CheckSteps"))
        {   
            prev = -1; 
            cout<<"\n Calling OutputPBandSE! at time"<<time;
            OutputPBandSE(time);
            check++;
        }
        prev++;
        
    }

    void UnsteadyAdvectionDiffusion::testComputeK(Array<OneD, Array<OneD, NekDouble> > &k)
    {

        // Number of fields (variables of the problem)
        int nVariables = m_fields.num_elements();
        
        //timestep
        double tchar = m_session->GetParameter("tchar");
       
        Array<OneD, NekDouble> temp = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
        Array<OneD, NekDouble> temp2 = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
        
        //pb[PBe2], 
        //N2starb*Pmaxb*Pba/s2s + N2starse*Pmaxse*Psea/s2s-e2btot //-z2btot not in Karin's code
        Vmath::Svtsvtp(nScaledDomainPts, N2starb*Pmaxb/s2s, PLPba, 1, N2starse*Pmaxse/s2s, PLPsea, 1, temp, 1);
        Vmath::Vsub(nScaledDomainPts, temp, 1, e2btot, 1, temp,1);
        Vmath::Vsub(nScaledDomainPts, temp, 1, z2btot, 1, temp,1);        
        //ke2on*e2*v2*temp
        Vmath::Vmul(nScaledDomainPts, temp, 1, fpInterp[FPe2], 1, temp, 1);
        Vmath::Smul(nScaledDomainPts, ke2on*v2, temp, 1, temp, 1);
        //-k2off*pb[PBe2] + kz2propbcat*pb[PBz2pro] // this second term is not in Karin's code
        Vmath::Svtsvtp(nScaledDomainPts, -k2off, pb[PBe2], 1, kz2propbcat, pb[PBz2pro], 1, temp2, 1);
        //Vmath::Smul(nScaledDomainPts, -ke2off, pb[PBe2], 1, temp2, 1);
        Vmath::Vadd(nScaledDomainPts, temp, 1, temp2, 1, temp, 1);
        //+ (kz5e2pbmin+kz5e2pbcat)*pb[PBz5e2]
        Vmath::Svtvp(nScaledDomainPts, (kz5e2pbmin+kz5e2pbcat), pb[PBz5e2], 1, temp, 1, temp, 1);
        //- kz5e2pbplu*pb[PBz5]*s5*pb[PBe2]
        Vmath::Vmul(nScaledDomainPts, pb[PBz5], 1, pb[PBe2], 1, temp2, 1);
        Vmath::Svtvp(nScaledDomainPts, -kz5e2pbplu*s5, temp2, 1, temp, 1, temp, 1);
        // + (kz8e2pbmin+kz8e2pbcat)*pb[PBz8e2]*s8/s2s
        Vmath::Svtvp(nScaledDomainPts, (kz8e2pbmin+kz8e2pbcat)*s8/s2s, pb[PBz8e2], 1, temp, 1, temp, 1);
        // - k15plu*PBz8*s8*pb[PBe2] 
        Vmath::Vmul(nScaledDomainPts, pb[PBz8], 1, pb[PBe2], 1, temp2, 1);
        Vmath::Svtvp(nScaledDomainPts,  -kz8e2pbplu*s8, temp2, 1, temp, 1, temp, 1);
        Vmath::Smul(nScaledDomainPts, tchar, temp, 1, k[3], 1);//>
        //cout<<"\n Inside testComputek  :: Max of k[PBe2] is :"<<Vmath::Vmax(nScaledDomainPts, &k[PBe2][0], 1)<<" min = "<<Vmath::Vmin(nScaledDomainPts, &k[PBe2][0], 1)<<" nvar is "<<3<< " PBe2 is"<< PBe2; 

    }

/**
 * @Brief Computes k1 and k2 for RK-2 stages for PB and SE reactions
 *
 * @param k  k1 or k2 depending on stage of RK-2
 */
    void UnsteadyAdvectionDiffusion::ComputeK(Array<OneD, Array<OneD, NekDouble> > &k)
    {
        int nVar = 0; // k has total nVar elements PB+SE
        
        // Number of fields (variables of the problem)
        int nVariables = m_fields.num_elements();
        
        //timestep
        double tchar = m_session->GetParameter("tchar");
       
        Array<OneD, NekDouble> temp = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
        Array<OneD, NekDouble> temp2 = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
        
        //Platelet-Bound species:
        
        //PBz5, 
        //N5b*Pmaxb*Pba/s5 + N5se*Pmaxse*Psea/s5-ze5btot        

        Vmath::Svtsvtp(nScaledDomainPts, N5b*Pmaxb/s5, PLPba, 1, N5se*Pmaxse/s5, PLPsea, 1, temp, 1);
        Vmath::Vsub(nScaledDomainPts, temp, 1, ze5btotInterp, 1, temp, 1);
    
        //k5on*FPINTERPz5*v5*temp
        Vmath::Vmul(nScaledDomainPts, temp, 1, fpInterp[FPz5], 1, temp, 1);
        Vmath::Smul(nScaledDomainPts, k5on*v5, temp, 1, temp, 1);
        
        //-k5off*pb[PBz5]
        Vmath::Svtvp(nScaledDomainPts, -k5off, pb[PBz5], 1, temp, 1, temp, 1);
        
        //+k5z5e10pbmin*pb[PBz5e10]*s10/s5 
        Vmath::Svtvp(nScaledDomainPts, kz5e10pbmin*s10/s5, pb[PBz5e10], 1, temp, 1, temp, 1);
      
        //-kz5e10pbplu*pb[PBz5]*PBe10*s10 
        Vmath::Vmul(nScaledDomainPts, pb[PBz5], 1, pb[PBe10], 1, temp2, 1);
        Vmath::Svtvp(nScaledDomainPts, -kz5e10pbplu*s10, temp2, 1, temp, 1, temp, 1);
       
        //+kz5e2pbmin*pb[PBz5]e2]*s2s/s5
        Vmath::Svtvp(nScaledDomainPts, kz5e2pbmin*s2s/s5, pb[PBz5e2], 1, temp, 1, temp, 1);
       
        //-kz5e2pbplu*pb[PBz5]*pb[PBe2]*s2s
        Vmath::Vmul(nScaledDomainPts, pb[PBz5], 1, pb[PBe2], 1, temp2, 1);
        Vmath::Svtvp(nScaledDomainPts, -kz5e2pbplu*s2s, temp2, 1, temp, 1, k[PBz5], 1);
        
        nVar++;
         
        // Galerkin project solution back to original space
        // m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, temp, pb[PBz5]);
        // for(int i = 0; i<nScaledDomainPts ; i++) 
        // {
        //     if(pb[PBz5][i]<0) pb[PBz5][i] = 0;
        // }

        //PBz8e2, 
        //kz8e2plu*PBz8*s2*pb[PBe2]
        Vmath::Smul(nScaledDomainPts, kz8e2pbplu*s2, pb[PBz8], 1, temp, 1);
        Vmath::Vmul(nScaledDomainPts, temp, 1, pb[PBe2], 1, temp, 1);
        //-(kz8e2min+kz8e2cat)*pb[PBz8e2]
        Vmath::Svtvp(nScaledDomainPts, -(kz8e2pbmin + kz8e2pbcat), pb[PBz8e2], 1, temp, 1, k[PBz8e2], 1);
        //Vmath::Smul(nScaledDomainPts, tchar, temp, 1, k[PBz8e2], 1);//>
        nVar++;
        // Galerkin project solution back to original space
        //m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, temp, pb[PBz8e2]);
        //for(int i = 0; i<nScaledDomainPts ; i++) 
        //{
        //    if(pb[PBz8e2][i]<0) pb[PBz8e2][i] = 0;
        //}

        //PBz5e2, 
        //kz5e2plu*pb[PBz5]*s5*pb[PBe2]
        Vmath::Smul(nScaledDomainPts, kz5e2pbplu*s5, pb[PBz5], 1, temp, 1);
        Vmath::Vmul(nScaledDomainPts, temp, 1, pb[PBe2], 1, temp, 1);
        //-(kz5e2min + kz5e2cat)pb[PBz5e2]
        Vmath::Svtvp(nScaledDomainPts, -(kz5e2pbmin + kz5e2pbcat), pb[PBz5e2], 1, temp, 1, k[PBz5e2], 1);
        //Vmath::Smul(nScaledDomainPts, tchar, temp, 1, k[PBz5e2], 1);//>
        nVar++;
        // Galerkin project solution back to original space
        //m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, temp, pb[PBz5e2]);
        //for(int i = 0; i<nScaledDomainPts ; i++) 
        //{
        //    if(pb[PBz5e2][i]<0) pb[PBz5e2][i] = 0;
        //}

        //PBe2, 
        //N2starb*Pmaxb*Pba/s2s + N2starse*Pmaxse*Psea/s2s-e2btot //-z2btot not in Karin's code
        Vmath::Svtsvtp(nScaledDomainPts, N2starb*Pmaxb/s2s, PLPba, 1, N2starse*Pmaxse/s2s, PLPsea, 1, temp, 1);
        
        // // Galerkin Projection e2btot back to non-scaled domain
        // m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, 
        //                                                 e2btotInterp, 
        //                                                 e2btot);        
        
        Vmath::Vsub(nScaledDomainPts, temp, 1, e2btotInterp, 1, temp,1);
        //Vmath::Vsub(nScaledDomainPts, temp, 1, z2btot, 1, temp,1);        
        //ke2on*e2*v2*temp
        Vmath::Vmul(nScaledDomainPts, temp, 1, fpInterp[FPe2], 1, temp, 1);
        Vmath::Smul(nScaledDomainPts, ke2on*v2, temp, 1, temp, 1);
        //-k2off*pb[PBe2] + kz2propbcat*pb[PBz2pro] // this second term is not in Karin's code
        //Vmath::Svtsvtp(nScaledDomainPts, -k2off, pb[PBe2], 1, kz2propbcat, pb[PBz2pro], 1, temp2, 1);
        Vmath::Smul(nScaledDomainPts, -ke2off, pb[PBe2], 1, temp2, 1);
        Vmath::Vadd(nScaledDomainPts, temp, 1, temp2, 1, temp, 1);
        //+ (kz5e2pbmin+kz5e2pbcat)*pb[PBz5e2]
        Vmath::Svtvp(nScaledDomainPts, (kz5e2pbmin+kz5e2pbcat), pb[PBz5e2], 1, temp, 1, temp, 1);
        //- kz5e2pbplu*pb[PBz5]*s5*pb[PBe2]
        Vmath::Vmul(nScaledDomainPts, pb[PBz5], 1, pb[PBe2], 1, temp2, 1);
        Vmath::Svtvp(nScaledDomainPts, -kz5e2pbplu*s5, temp2, 1, temp, 1, temp, 1);
        // + (kz8e2pbmin+kz8e2pbcat)*pb[PBz8e2]*s8/s2s
        Vmath::Svtvp(nScaledDomainPts, (kz8e2pbmin+kz8e2pbcat)*s8/s2s, pb[PBz8e2], 1, temp, 1, temp, 1);
        // - k15plu*PBz8*s8*pb[PBe2] 
        Vmath::Vmul(nScaledDomainPts, pb[PBz8], 1, pb[PBe2], 1, temp2, 1);
        Vmath::Svtvp(nScaledDomainPts,  -kz8e2pbplu*s8, temp2, 1, temp, 1, k[PBe2], 1);
        //Vmath::Smul(nScaledDomainPts, tchar, temp, 1, k[PBe2], 1);//>
        nVar++;
       // Galerkin project solution back to original space
        // m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, temp, pb[PBe2]);
        // for(int i = 0; i<nScaledDomainPts ; i++) 
        // {
        //     if(k[PBe2][i]<0) k[PBe2][i] = 0;
        // }
   
        //PBz5e10, 
        //kz5e2pbplu*pb[PBz5]*s5*PBe10
        Vmath::Vmul(nScaledDomainPts, pb[PBz5], 1, pb[PBe10], 1, temp, 1);
        Vmath::Smul(nScaledDomainPts, kz5e2pbplu*s5, temp, 1, temp, 1);
        //-(kz5e2pbmin+kz5e2pbcat)*pb[PBz5e10]
        Vmath::Svtvp(nScaledDomainPts, -(kz5e2pbmin+kz5e2pbcat), pb[PBz5e10], 1, temp, 1, k[PBz5e10], 1);
//        Vmath::Smul(nScaledDomainPts, tchar, temp, 1, k[PBz5e10], 1);//>
        nVar++;
//        cout<<"\n 11. Inside Computek  :: Max of k[PBe2] is :"<<Vmath::Vmax(nScaledDomainPts, k[PBe2], 1)<<" min = "<<Vmath::Vmin(nScaledDomainPts, k[PBe2], 1)<<" nvar is "<<nVar<< " PBe2 is"<< PBe2; 

        // Galerkin project solution back to original space
        // m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, temp, pb[PBz5e10]);
        // for(int i = 0; i<nScaledDomainPts ; i++) 
        // {
        //     if(pb[PBz5e10][i]<0) pb[PBz5e10][i] = 0;
        // }

    
        //PBe5, 
        //N5b*Pmaxb*Pba/s5 + N5se*Pmaxse*Psea/s5-ze5btot
        Vmath::Svtsvtp(nScaledDomainPts, N5b*Pmaxb/s5, PLPba, 1, N5se*Pmaxse/s5, PLPsea, 1, temp, 1);
        Vmath::Vsub(nScaledDomainPts, temp, 1, ze5btotInterp, 1, temp, 1);
        //k5on*FPe5*v5*temp
        Vmath::Vmul(nScaledDomainPts, temp, 1, fpInterp[FPe5], 1, temp, 1);
        Vmath::Smul(nScaledDomainPts, k5on*v5, temp, 1, temp, 1);
        //-k5off*pb[PBe5]
        Vmath::Svtvp(nScaledDomainPts, -k5off, pb[PBe5], 1, temp, 1, temp, 1);
        //+kz5e10pbcat*pb[PBz5e10]*s10/s5 
        Vmath::Svtvp(nScaledDomainPts, kz5e10pbcat*s10/s5, pb[PBz5e10], 1, temp, 1, temp, 1);
        //+kz5e2pbcat*pb[PBz5]e2]*s2s/s5 
        Vmath::Svtvp(nScaledDomainPts, kz5e2pbcat*s2s/s5, pb[PBz5e2], 1, temp, 1, temp, 1);
        //+kpromin*pb[PBpro]
        Vmath::Svtvp(nScaledDomainPts, kpromin, pb[PBpro], 1, temp, 1, temp, 1);
        //-kproplu*pb[PBe5]*PBe10*s10
        Vmath::Vmul(nScaledDomainPts, pb[PBe5], 1, pb[PBe10], 1, temp2, 1);
        Vmath::Svtvp(nScaledDomainPts, -kproplu*s10, temp2, 1, temp, 1, temp, 1);
        //+kapce5pbmin*PBapce5
        Vmath::Svtvp(nScaledDomainPts, kapce5pbmin, pb[PBapce5], 1, temp, 1, temp, 1);
        // -kapce5pbplu*FPapc*v2*pb[PBe5]
        Vmath::Vmul(nScaledDomainPts, pb[PBe5], 1, fpInterp[FPapc], 1, temp2, 1);
        Vmath::Svtvp(nScaledDomainPts, -kapce5pbplu*v2, temp2, 1, temp, 1, k[PBe5], 1);
        //Vmath::Smul(nScaledDomainPts, tchar, temp, 1, k[PBe5], 1);//>
        nVar++;
        // Galerkin project solution back to original space
        // m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, temp, pb[PBe5]);
        // for(int i = 0; i<nScaledDomainPts ; i++) 
        // {
        //     if(pb[PBe5][i]<0) pb[PBe5][i] = 0;
        // }

    
        //PBapce5, 
        //kapce5pbplu*FPapc*v2*pb[PBe5]
        Vmath::Vmul(nScaledDomainPts, fpInterp[FPapc], 1, pb[PBe5], 1, temp, 1);
        Vmath::Smul(nScaledDomainPts, kapce5pbplu*v2, temp, 1, temp, 1);
        //-(kapce5pbmin+kapce5pbcat)*PBapce5
        Vmath::Svtvp(nScaledDomainPts, -(kapce5pbmin+kapce5pbcat), pb[PBapce5], 1, temp, 1, k[PBapce5], 1);
//        Vmath::Smul(nScaledDomainPts, tchar, temp, 1, k[PBapce5], 1);//>
        nVar++;
        // Galerkin project solution back to original space
        // m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, temp, pb[PBapce5]);
        // for(int i = 0; i<nScaledDomainPts ; i++) 
        // {
        //     if(pb[PBapce5][i]<0) pb[PBapce5][i] = 0;
        // }

    
        //pb[PBpro],
        //kproplu*pb[PBe5]*PBe10*s10
        Vmath::Vmul(nScaledDomainPts, pb[PBe5], 1, pb[PBe10], 1, temp, 1);
        Vmath::Smul(nScaledDomainPts, kproplu*s10, temp, 1, temp, 1);
        //-kpromin*pb[PBpro]
        Vmath::Svtvp(nScaledDomainPts, -kpromin, pb[PBpro], 1, temp, 1, temp, 1);
        //+(kz2propbmin+kz2propbcat)*pb[PBz2pro] 
        Vmath::Svtvp(nScaledDomainPts, (kz2propbmin+kz2propbcat), pb[PBz2pro], 1, temp, 1, temp, 1);
        //-kz2propbplu*PBz2*s2*pb[PBpro]
        Vmath::Vmul(nScaledDomainPts, pb[PBz2], 1, pb[PBpro], 1, temp2, 1);
        Vmath::Svtvp(nScaledDomainPts, -kz2propbplu*s2, temp2, 1, temp, 1, k[PBpro], 1);
        //Vmath::Smul(nScaledDomainPts, tchar, temp, 1, k[PBpro], 1);//>
        nVar++;
        // Galerkin project solution back to original space
        // m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, temp, pb[PBpro]);
        // for(int i = 0; i<nScaledDomainPts ; i++) 
        // {
        //     if(pb[PBpro][i]<0) pb[PBpro][i] = 0;
        // }
    
        //pb[PBz2pro], 
        //kz2propbplu*PBz2*s2*pro
        Vmath::Vmul(nScaledDomainPts, pb[PBz2], 1, pb[PBpro], 1, temp, 1);
        Vmath::Smul(nScaledDomainPts, kz2propbplu*s2, temp, 1, temp, 1);
        //-(kz2propbmin+kz2propbcat)*pb[PBz2pro]
        Vmath::Svtvp(nScaledDomainPts, -(kz2propbmin+kz2propbcat), pb[PBz2pro], 1, temp, 1, k[PBz2pro], 1);//>
//        Vmath::Smul(nScaledDomainPts, tchar, temp, 1, k[PBz2pro], 1);//>
        nVar++;
        // Galerkin project solution back to original space
        // m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, temp, pb[PBz2pro]);
        // for(int i = 0; i<nScaledDomainPts ; i++) 
        // {
        //     if(pb[PBz2pro][i]<0) pb[PBz2pro][i] = 0;
        // }

        //pb[PBz10tenstar], 
        //kz10tenpbplu*PBz10*s10*tenstar
        Vmath::Vmul(nScaledDomainPts, pb[PBz10], 1, pb[PBtenstar], 1, temp, 1);
        Vmath::Smul(nScaledDomainPts, kz10tenpbplu*s10, temp, 1, temp, 1);
        //-(kz1tenpbmin+kz10tenpbcat)*pb[PBz10tenstar]
        Vmath::Svtvp(nScaledDomainPts, -(kz10tenpbmin+kz10tenpbcat), pb[PBz10tenstar], 1, temp, 1, k[PBz10tenstar], 1);
//        Vmath::Smul(nScaledDomainPts, tchar, temp, 1, k[PBz10tenstar], 1);//>
        nVar++;
        // Galerkin project solution back to original space
        // m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, temp, pb[PBz10tenstar]);   
        // for(int i = 0; i<nScaledDomainPts ; i++) 
        // {
        //     if(pb[PBz10tenstar][i]<0) pb[PBz10tenstar][i] = 0;
        // }
    
        //PBz8e10, 
        //kz8e10pbplu*PBz8*PBe10*s10
        Vmath::Vmul(nScaledDomainPts, pb[PBz8], 1, pb[PBe10], 1, temp, 1);
        Vmath::Smul(nScaledDomainPts, kz8e10pbplu*s10, temp, 1, temp, 1);
        //-(kz8e10pbmin+kz8e10pbcat)*PBz8e10
        Vmath::Svtvp(nScaledDomainPts, -(kz8e10pbmin+kz8e10pbcat), pb[PBz8e10], 1, temp, 1, k[PBz8e10], 1);
        //Vmath::Smul(nScaledDomainPts, tchar, temp, 1, k[PBz8e10], 1);//>
        nVar++;
        // Galerkin project solution back to original space
        // m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, temp, pb[PBz8e10]);   
        // for(int i = 0; i<nScaledDomainPts ; i++) 
        // {
        //     if(pb[PBz8e10][i]<0) pb[PBz8e10][i] = 0;
        // }

        //PBz10ten, 
        //kz10tenpbplu*PBz10*s10*PBten
        Vmath::Vmul(nScaledDomainPts, pb[PBz10], 1, pb[PBten], 1, temp, 1);
        Vmath::Smul(nScaledDomainPts, kz10tenpbplu*s10, temp, 1, temp, 1);
        //-(kz10tenpbmin+kz10tenpbcat)*PBz10ten
        Vmath::Svtvp(nScaledDomainPts, -(kz10tenpbmin+kz10tenpbcat), pb[PBz10ten], 1, temp, 1, k[PBz10ten], 1);
//        Vmath::Smul(nScaledDomainPts, tchar, temp, 1, k[PBz10ten], 1);//>
        nVar++;
        
        // Galerkin project solution back to original space
        // m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, temp, pb[PBz10ten]);   
        // for(int i = 0; i<nScaledDomainPts ; i++) 
        // {
        //     if(pb[PBz10ten][i]<0) pb[PBz10ten][i] = 0;
        // }
    
        //PBe10, 
        //N10b*Pmaxb*Pba/s10 + N10se*Pmaxse*Psea/s10
        Vmath::Svtsvtp(nScaledDomainPts, N10b*Pmaxb/s10, PLPba, 1, N10se*Pmaxse/s10, PLPsea, 1, temp, 1);
       
        //-ze10btot
    
        // // Galerkin Projection ze10btot back to non-scaled domain
        // m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, 
        //                                                 ze10btotInterp, 
        //                                                 ze10btot);
                
        Vmath::Vsub(nScaledDomainPts, temp, 1, ze10btotInterp, 1, temp, 1);
       
        //k10on*FPe10*v10*temp
        Vmath::Vmul(nScaledDomainPts, temp, 1, fpInterp[FPe10], 1, temp, 1);
        //cout<<" \n temp3 is : "<<Vmath::Vmax(nScaledDomainPts, temp, 1)<<" min = "<<Vmath::Vmin(nScaledDomainPts, temp, 1);

        Vmath::Smul(nScaledDomainPts, k10on*v10, temp, 1, temp, 1);
      
        //-k10off*PBe10
        Vmath::Svtvp(nScaledDomainPts, -k10off, pb[PBe10], 1, temp, 1, temp, 1);
      
        //+kz10tenpbcat*PBz10ten*s8/s10  
        Vmath::Svtvp(nScaledDomainPts, kz10tenpbcat*s8/s10, pb[PBz10ten], 1, temp, 1, temp, 1);
     
        //+(kz5e10pbmin+kz5e10pbcat)*pb[PBz5e10]
        Vmath::Svtvp(nScaledDomainPts, kz5e10pbmin + kz5e10pbcat, pb[PBz5e10], 1, temp, 1, temp, 1);
     
        //+(kz8e10pbmin+kz8e10pbcat)*PBz8e10*s8/s10 
        Vmath::Svtvp(nScaledDomainPts, (kz8e10pbmin + kz8e10pbcat)*s8/s10 , pb[PBz8e10], 1, temp, 1, temp, 1);
    
        //-kz5e10pbplu*pb[PBz5]*s5*PBe10
        Vmath::Vmul(nScaledDomainPts, pb[PBz5], 1, pb[PBe10], 1, temp2, 1);
    
        Vmath::Svtvp(nScaledDomainPts, -kz5e10pbplu*s5, temp2, 1, temp, 1, temp, 1);
    
        //-kz8e10pbplu*PBz8*s8*PBe10
        Vmath::Vmul(nScaledDomainPts, pb[PBz8], 1, pb[PBe10], 1, temp2, 1);

        Vmath::Svtvp(nScaledDomainPts, -kz8e10pbplu*s8, temp2, 1, temp, 1, temp, 1);

        //+kpromin*pb[PBpro]*s5/s10 
        Vmath::Svtvp(nScaledDomainPts, kpromin*s5/s10, pb[PBpro], 1, temp, 1, temp, 1);

        //-kproplu*pb[PBe5]*s5*PBe10
        Vmath::Vmul(nScaledDomainPts, pb[PBe5], 1, pb[PBe10], 1, temp2, 1);

        Vmath::Svtvp(nScaledDomainPts, -kproplu*s5, temp2, 1, temp, 1, temp, 1);

        //+kz10tenpbcat*pb[PBz10tenstar]*s8/s10
        Vmath::Svtvp(nScaledDomainPts, kz10tenpbcat*s8/s10, pb[PBz10tenstar], 1, temp, 1, k[PBe10], 1);

        //cout<<" \n temp16 is : "<<Vmath::Vmax(nScaledDomainPts, temp, 1)<<" min = "<<Vmath::Vmin(nScaledDomainPts, temp, 1);
//        Vmath::Smul(nScaledDomainPts, tchar, temp, 1, k[PBe10], 1);//>
         //cout<<"\n Max of k[PBe10] is :"<<Vmath::Vmax(nScaledDomainPts, k[PBe10], 1)<<" min = "<<Vmath::Vmin(nScaledDomainPts, k[PBe10], 1);
        //cout<<" Max of pb[PBe10] is :"<<Vmath::Vmax(nScaledDomainPts, pb[PBe10], 1)<<" min = "<<Vmath::Vmin(nScaledDomainPts, pb[PBe10], 1);
        nVar++;
        
        // Galerkin project solution back to original space
        // m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, temp, pb[PBe10]);   
        // for(int i = 0; i<nScaledDomainPts ; i++) 
        // {
        //     if(pb[PBe10][i]<0) pb[PBe10][i] = 0;
        // }
    
        //PBz10, 
        //N10b*Pmaxb*Pba/s10 + N10se*Pmaxse*Psea/s10
        Vmath::Svtsvtp(nScaledDomainPts, N10b*Pmaxb/s10, PLPba, 1, N10se*Pmaxse/s10, PLPsea, 1, temp, 1);
        //-ze10btot
        Vmath::Vsub(nScaledDomainPts, temp, 1, ze10btotInterp, 1, temp, 1);
        //k10on*FPz10*v10*temp
        Vmath::Vmul(nScaledDomainPts, temp, 1, fpInterp[FPz10], 1, temp, 1);
        Vmath::Smul(nScaledDomainPts, k10on*v10, temp, 1, temp, 1);
        //-k10zoff*PBz10
        Vmath::Svtvp(nScaledDomainPts, k10zoff, pb[PBz10], 1, temp, 1, temp, 1);
        //+kz10tenpbmin*PBz10ten*s8/s10 
        Vmath::Svtvp(nScaledDomainPts, kz10tenpbmin*s8/s10, pb[PBz10ten], 1, temp, 1, temp, 1);
        //-kz10tenpbplu*PBz10*PBten*s8 
        Vmath::Vmul(nScaledDomainPts, pb[PBz10], 1, pb[PBten], 1, temp2, 1);
        Vmath::Svtvp(nScaledDomainPts, -kz10tenpbplu*s8, temp2, 1, temp, 1, temp, 1);
        //+kz10tenpbmin*pb[PBz10tenstar]*s8/s10
        Vmath::Svtvp(nScaledDomainPts, kz10tenpbmin*s8/s10, pb[PBz10tenstar], 1, temp, 1, temp, 1);
        //-kz10tenpbplu*PBz10*PBtenstar*s8
        Vmath::Vmul(nScaledDomainPts, pb[PBz10], 1, pb[PBtenstar], 1, temp2, 1);
        Vmath::Svtvp(nScaledDomainPts, -kz10tenpbplu*s8, temp2, 1, temp, 1, k[PBz10], 1);
//        Vmath::Smul(nScaledDomainPts, tchar, temp, 1, k[PBz10], 1);//>
        nVar++;
        // Galerkin project solution back to original space
        // m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, temp, pb[PBz10]);   
        // for(int i = 0; i<nScaledDomainPts ; i++) 
        // {
        //     if(pb[PBz10][i]<0) pb[PBz10][i] = 0;
        // }
    
        //PBe8, 
        //N8b*Pmaxb*Pba/s8 + N8se*Pmaxse*Psea/s8
        Vmath::Svtsvtp(nScaledDomainPts, N8b*Pmaxb/s8, PLPba, 1, N8se*Pmaxse/s8, PLPsea, 1, temp, 1);
        //-ze8btot
        
        // // Galerkin Projection ze8btot back to non-scaled domain
        // m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, 
        //                                                 ze8btotInterp, 
        //                                                 ze8btot);       
        
        Vmath::Vsub(nScaledDomainPts, temp, 1, ze8btotInterp, 1, temp, 1);  
        //k8on*FPe8*v8*temp
        Vmath::Vmul(nScaledDomainPts, fpInterp[FPe8], 1, temp, 1, temp, 1);
        Vmath::Smul(nScaledDomainPts, k8on*v8, temp, 1, temp, 1);
        //-k8off*e8ba 
        Vmath::Svtvp(nScaledDomainPts, -k8off, pb[PBe8], 1, temp, 1, temp, 1);
        //+kz8e10pbcat*PBz8e10 
        Vmath::Svtvp(nScaledDomainPts, kz8e10pbcat, pb[PBz8e10], 1, temp, 1, temp, 1);
        //+kz8e2pbcat*pb[PBz8e2]
        Vmath::Svtvp(nScaledDomainPts, kz8e2pbcat, pb[PBz8e2], 1, temp, 1, temp, 1);
        //+ktenmin*PBten 
        Vmath::Svtvp(nScaledDomainPts, ktenmin, pb[PBten], 1, temp, 1, temp, 1);
        //-ktenplu*PBe8*PBe9*s9 
        Vmath::Vmul(nScaledDomainPts, pb[PBe8], 1, pb[PBe9], 1, temp2, 1);
        Vmath::Svtvp(nScaledDomainPts, -ktenplu*s9, temp2, 1, temp, 1, temp, 1);
        //+kapce8pbmin*PBapce8 
        Vmath::Svtvp(nScaledDomainPts, kapce8pbmin, pb[PBapce8], 1, temp, 1, temp, 1);
        //-kapce8pbplu*FPapc*v2*PBe8
        Vmath::Vmul(nScaledDomainPts, pb[PBe8], 1, fpInterp[FPapc], 1, temp2, 1);
        Vmath::Svtvp(nScaledDomainPts, -kapce8pbplu*v2, temp2, 1, temp, 1, temp, 1);
        //-k10plu*PBe8*PBe9star*s91 
        Vmath::Vmul(nScaledDomainPts, pb[PBe8], 1, pb[PBe9star], 1, temp2, 1);
        Vmath::Svtvp(nScaledDomainPts, -kapce8pbplu*s91, temp2, 1, temp, 1, temp, 1);
        //+ktenmin*PBtenstar
        Vmath::Svtvp(nScaledDomainPts, ktenmin, pb[PBtenstar], 1, temp, 1, k[PBe8], 1);
//        Vmath::Smul(nScaledDomainPts, tchar, temp, 1, k[PBe8], 1);//>
        nVar++;
        // Galerkin project solution back to original space
        // m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, temp, pb[PBe8]);   
        // for(int i = 0; i<nScaledDomainPts ; i++) 
        // {
        //     if(pb[PBe8][i]<0) pb[PBe8][i] = 0;
        // }
    
        //PBtenstar, 
        //ktenplu*PBe8*PBe9star*s91
        Vmath::Vmul(nScaledDomainPts, pb[PBe8], 1, pb[PBe9star], 1, temp, 1);
        Vmath::Smul(nScaledDomainPts, ktenplu*s91, temp, 1, temp, 1);
        //-ktenmin*PBtenstar 
        Vmath::Svtvp(nScaledDomainPts, -ktenmin, pb[PBtenstar], 1, temp, 1, temp, 1);
        //+(kz10tenpbmin+kz10tenpbcat)*pb[PBz10tenstar] 
        Vmath::Svtvp(nScaledDomainPts, (kz10tenpbmin+kz10tenpbcat), pb[PBtenstar], 1, temp, 1, temp, 1);
        //-kz10tenpbplu*PBz10*s10*PBtenstar
        Vmath::Vmul(nScaledDomainPts, pb[PBz10], 1, pb[PBtenstar], 1, temp2, 1);
        Vmath::Svtvp(nScaledDomainPts, -kz10tenpbplu*s10, temp2, 1, temp, 1, k[PBtenstar], 1);
//        Vmath::Smul(nScaledDomainPts, tchar, temp, 1, k[PBtenstar], 1);//>
        nVar++;
        // Galerkin project solution back to original space
        // m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, temp, pb[PBtenstar]); 
        // for(int i = 0; i<nScaledDomainPts ; i++) 
        // {
        //     if(pb[PBtenstar][i]<0) pb[PBtenstar][i] = 0;
        // }   
    
        //PBten,
        //ktenplu*PBe8*PBe9*s9 
        Vmath::Vmul(nScaledDomainPts, pb[PBe8], 1, pb[PBe9], 1, temp, 1);
        Vmath::Smul(nScaledDomainPts, ktenplu*s9, temp, 1, temp, 1);
        //-ktenmin*ten 
        Vmath::Svtvp(nScaledDomainPts, -ktenmin, pb[PBten], 1, temp, 1, temp, 1);
        //+(kz10tenpbmin+kz10tenpbcat)*PBz10ten 
        Vmath::Svtvp(nScaledDomainPts, (kz10tenpbmin+kz10tenpbcat), pb[PBz10ten], 1, temp, 1, temp, 1);
        //-kz10tenpbplu*PBz10*s10*PBten 
        Vmath::Vmul(nScaledDomainPts, pb[PBz10], 1, pb[PBten], 1, temp2, 1);
        Vmath::Svtvp(nScaledDomainPts, -kz10tenpbplu*s10, temp2, 1, temp, 1, k[PBten], 1);
//        Vmath::Smul(nScaledDomainPts, tchar, temp, 1, k[PBten], 1);//>
        nVar++;
        // Galerkin project solution back to original space
        // m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, temp, pb[PBten]);    
        // for(int i = 0; i<nScaledDomainPts ; i++) 
        // {
        //     if(pb[PBten][i]<0) pb[PBten][i] = 0;
        // }
    
        //PBz8, 
        //N8b*Pmaxb*Pba/s8 + N8se*Pmaxse*Psea/s8
        Vmath::Svtsvtp(nScaledDomainPts, N8b*Pmaxb/s8, PLPba, 1, N8se*Pmaxse/s8, PLPsea, 1, temp, 1);
        // cout<<" \n temp1 is : "<<Vmath::Vmax(nScaledDomainPts, temp, 1)<<" min = "<<Vmath::Vmin(nScaledDomainPts, temp, 1);

        //-ze8btot
        Vmath::Vsub(nScaledDomainPts, temp, 1, ze8btotInterp, 1, temp, 1);  

        // cout<<" \n temp2 is : "<<Vmath::Vmax(nScaledDomainPts, temp, 1)<<" min = "<<Vmath::Vmin(nScaledDomainPts, temp, 1);

        //k8on*z8*v8*temp
        Vmath::Vmul(nScaledDomainPts, fpInterp[FPz8], 1, temp, 1, temp, 1);
        
        // cout<<" \n temp3 is : "<<Vmath::Vmax(nScaledDomainPts, temp, 1)<<" min = "<<Vmath::Vmin(nScaledDomainPts, temp, 1);

        Vmath::Smul(nScaledDomainPts, k8on*v8, temp, 1, temp, 1);
        
        // cout<<" \n temp4 is : "<<Vmath::Vmax(nScaledDomainPts, temp, 1)<<" min = "<<Vmath::Vmin(nScaledDomainPts, temp, 1);

        //-k8off*PBz8
        Vmath::Svtvp(nScaledDomainPts, -k8off, pb[PBz8], 1, temp, 1, temp, 1);
        
        // cout<<" \n temp5 is : "<<Vmath::Vmax(nScaledDomainPts, temp, 1)<<" min = "<<Vmath::Vmin(nScaledDomainPts, temp, 1);

        //+kz8e10pbmin*PBz8e10
        Vmath::Svtvp(nScaledDomainPts, -kz8e10pbmin, pb[PBz8e10], 1, temp, 1, temp, 1);
        // cout<<" \n temp6 is : "<<Vmath::Vmax(nScaledDomainPts, temp, 1)<<" min = "<<Vmath::Vmin(nScaledDomainPts, temp, 1);
        //-kz8e10pbplu*PBz8*PBe10*s10 
        Vmath::Vmul(nScaledDomainPts, pb[PBz8], 1, pb[PBe10], 1, temp2, 1);
        
        // cout<<" \n temp7 is : "<<Vmath::Vmax(nScaledDomainPts, temp2, 1)<<" min = "<<Vmath::Vmin(nScaledDomainPts, temp2, 1);

        Vmath::Svtvp(nScaledDomainPts, -kz8e10pbplu*s10, temp2, 1, temp, 1, temp, 1);
        
        // cout<<" \n temp8 is : "<<Vmath::Vmax(nScaledDomainPts, temp, 1)<<" min = "<<Vmath::Vmin(nScaledDomainPts, temp, 1);

        //+kz8e2pbmin*pb[PBz8e2] 
        Vmath::Svtvp(nScaledDomainPts, kz8e2pbmin, pb[PBz8e2], 1, temp, 1, temp, 1);
        
        // cout<<" \n temp9 is : "<<Vmath::Vmax(nScaledDomainPts, temp, 1)<<" min = "<<Vmath::Vmin(nScaledDomainPts, temp, 1);

        //-kz8e2pbplu*PBz8*pb[PBe2]*s2s
        Vmath::Vmul(nScaledDomainPts, pb[PBz8], 1, pb[PBe2], 1, temp2, 1);
        
        // cout<<" \n temp10 is : "<<Vmath::Vmax(nScaledDomainPts, temp2, 1)<<" min = "<<Vmath::Vmin(nScaledDomainPts, temp2, 1);

        Vmath::Svtvp(nScaledDomainPts, -kz8e2pbplu*s2s, temp2, 1, temp, 1,  k[PBz8], 1);
        
        // cout<<" \n temp11 is : "<<Vmath::Vmax(nScaledDomainPts, temp, 1)<<" min = "<<Vmath::Vmin(nScaledDomainPts, temp, 1);

//        Vmath::Smul(nScaledDomainPts, tchar, temp, 1, k[PBz8], 1);//>
        // cout<<"\n 1. Inside Computek  :: Max of pb[PBz8] is :"<<Vmath::Vmax(nScaledDomainPts, pb[PBz8], 1)<<" min = "<<Vmath::Vmin(nScaledDomainPts, pb[PBz8], 1); 
        //cout<<" Max of PLPba is :"<<Vmath::Vmax(nScaledDomainPts, PLPba, 1)<<" min = "<<Vmath::Vmin(nScaledDomainPts, PLPba, 1);
        //cout<<" Max of PLPsea is :"<<Vmath::Vmax(nScaledDomainPts, PLPsea, 1)<<" min = "<<Vmath::Vmin(nScaledDomainPts, PLPsea, 1);
       // cout<<" Max of fpInterp[FPz8] is :"<<Vmath::Vmax(nScaledDomainPts, fpInterp[FPz8], 1)<<" min = "<<Vmath::Vmin(nScaledDomainPts, fpInterp[FPz8], 1);
        //cout<<" Max of pb[PBz8e10] is :"<<Vmath::Vmax(nScaledDomainPts, pb[PBz8e10], 1)<<" min = "<<Vmath::Vmin(nScaledDomainPts, pb[PBz8e10], 1);
        //cout<<" Max of pb[PBz8e2] is :"<<Vmath::Vmax(nScaledDomainPts, pb[PBz8e2], 1)<<" min = "<<Vmath::Vmin(nScaledDomainPts, pb[PBz8e2], 1);
        // cout<<" Max of pb[PBe10] is :"<<Vmath::Vmax(nScaledDomainPts, pb[PBe10], 1)<<" min = "<<Vmath::Vmin(nScaledDomainPts, pb[PBe10], 1);
        
        // cout<<" Max of k[PBz8] is :"<<Vmath::Vmax(nScaledDomainPts, k[PBz8], 1)<<" min = "<<Vmath::Vmin(nScaledDomainPts, k[PBz8], 1);
        nVar++;
        // Galerkin project solution back to original space
        // m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, temp, pb[PBz8]);    
        // for(int i = 0; i<nScaledDomainPts ; i++) 
        // {
        //     if(pb[PBz8][i]<0) pb[PBz8][i] = 0;
        // }
    
        //PBapce8
        //kapce8pbplu*FPapc*v2*PBe8
        Vmath::Vmul(nScaledDomainPts, fpInterp[FPapc], 1, pb[PBe8], 1, temp, 1);
        Vmath::Smul(nScaledDomainPts, kapce8pbplu*v2, temp, 1, temp, 1);
        //-(kapce8pbmin+kapce8pbcat)*PBapce8
        Vmath::Svtvp(nScaledDomainPts, -(kapce8pbmin+kapce8pbcat), pb[PBapce8], 1, temp, 1, k[PBapce8], 1);
        //Vmath::Smul(nScaledDomainPts, tchar, temp, 1, k[PBapce8], 1);//>
        nVar++;
        // Galerkin project solution back to original space
        // m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, temp, pb[PBapce8]);    
        // for(int i = 0; i<nScaledDomainPts ; i++) 
        // {
        //     if(pb[PBapce8][i]<0) pb[PBapce8][i] = 0;
        // }
  
        //PBz2
        //k2on*z2*v2  
        Vmath::Smul(nScaledDomainPts, k2on*v2, fpInterp[FPz2], 1, temp2, 1);
        //(N2b*Pmaxb*Pba/s2 + N2se*Pmaxse*Psea/s2 - z2btot // -e2btot ) last term not in karin's code
        Vmath::Svtsvtp(nScaledDomainPts, N2b*Pmaxb/s2, PLPba, 1, N2se*Pmaxse/s2, PLPsea, 1, temp, 1);

        // // Galerkin Projection z2btot back to non-scaled domain
        // m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, 
        //                                                 z2btotInterp, 
        //                                                 z2btot);    
                
        Vmath::Vsub(nScaledDomainPts, temp, 1, z2btotInterp, 1, temp, 1);
        Vmath::Vmul(nScaledDomainPts, temp, 1, temp2, 1, temp, 1);
        //-k2off*PBz2 + k7min*pb[PBz2pro]*s5/s2 
        Vmath::Svtsvtp(nScaledDomainPts, -k2off, pb[PBz2], 1, kz2propbmin*s5/s2, pb[PBz2pro], 1, temp2, 1);
        Vmath::Vadd(nScaledDomainPts, temp, 1, temp2, 1, temp, 1);
        // - k7plu*PBz2*pb[PBpro]*s5
        Vmath::Vmul(nScaledDomainPts, pb[PBz2], 1, pb[PBpro], 1, temp2, 1);
        Vmath::Svtvp(nScaledDomainPts, -kz2propbplu*s5, temp2, 1, temp, 1,  k[PBz2], 1);
//        Vmath::Smul(nScaledDomainPts, tchar, temp, 1, k[PBz2], 1);//>
        nVar++;
        // Galerkin project solution back to original space
        // m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, temp, pb[PBz2]);
        // for(int i = 0; i<nScaledDomainPts ; i++) 
        // {
        //     if(pb[PBz2][i]<0) pb[PBz2][i] = 0;
        // }

        //PBe9
        //N9b*Pmaxb*Pba/s9 + N9se*Pmaxse*Psea/s9
        Vmath::Svtsvtp(nScaledDomainPts, N9b*Pmaxb/s9, PLPba, 1, N9se*Pmaxse/s9, PLPsea, 1, temp, 1);
        //-ze9btot
        // // Galerkin Projection ze5btot back to non-scaled domain
        // m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, 
        //                                                 ze9btotInterp, 
        //                                                 ze9btot);    
                
        Vmath::Vsub(nScaledDomainPts, temp, 1, ze9btotInterp, 1, temp, 1);  
        //k9on*FPe9*v9*temp
        Vmath::Vmul(nScaledDomainPts, temp, 1, fpInterp[FPe9], 1, temp, 1);
        Vmath::Smul(nScaledDomainPts, k9on*v9, temp, 1, temp, 1);
        //-k9off*PBe9
        Vmath::Svtvp(nScaledDomainPts, -k9off, pb[PBe9], 1, temp, 1, temp, 1);   
        //+ktenmin*PBten*s8/s9 
        Vmath::Svtvp(nScaledDomainPts, ktenmin*s8/s9, pb[PBten], 1, temp, 1, temp, 1);
        //-ktenplu*PBe8*s8*PBe9 
        Vmath::Vmul(nScaledDomainPts, pb[PBe8], 1, pb[PBe9], 1, temp2, 1);
        Vmath::Svtvp(nScaledDomainPts, -ktenplu*s8, temp2, 1, temp, 1,  k[PBe9], 1);
//        Vmath::Smul(nScaledDomainPts, tchar, temp, 1, k[PBe9], 1);//>
       
        nVar++;
        // Galerkin project solution back to original space
        // m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, temp, pb[PBe9]);
        // for(int i = 0; i<nScaledDomainPts ; i++) 
        // {
        //     if(pb[PBe9][i]<0) pb[PBe9][i] = 0;
        // }

        //PBe9star
        //N9starb*Pmaxb*Pba + N9starse*Pmaxse*Psea - PBe9star*s91 - PBtenstar*s8 - pb[PBz10tenstar]*s8
        Vmath::Svtsvtp(nScaledDomainPts, N9starb*Pmaxb, PLPba, 1, N9starse*Pmaxse, PLPsea, 1, temp, 1);
        Vmath::Svtsvtp(nScaledDomainPts, -s91, pb[PBe9star], 1, -s8, pb[PBtenstar], 1, temp2, 1);
        Vmath::Vadd(nScaledDomainPts, temp, 1, temp2, 1, temp, 1);
        Vmath::Svtvp(nScaledDomainPts, -s8, pb[PBz10tenstar], 1, temp, 1, temp, 1);
        //k9on*FPe9*v9/s91*temp
        Vmath::Vmul(nScaledDomainPts, fpInterp[FPe9], 1, temp, 1, temp, 1);
        Vmath::Smul(nScaledDomainPts, k9on*v9/s91, temp, 1, temp, 1);
        //-k9off*PBe9star
        Vmath::Svtvp(nScaledDomainPts, -k9off, pb[PBe9star], 1, temp, 1, temp, 1);
        //+ktenmin*PBtenstar*s8/s91 
        Vmath::Svtvp(nScaledDomainPts, ktenmin*s8/s91, pb[PBtenstar], 1, temp, 1, temp, 1);
        //-ktenplu*PBe8*s8*PBe9star
        Vmath::Vmul(nScaledDomainPts, pb[PBe8], 1, pb[PBe9star], 1, temp2, 1);
        Vmath::Svtvp(nScaledDomainPts, -ktenplu*s8, temp2, 1, temp, 1,  k[PBe9star], 1);
//        Vmath::Smul(nScaledDomainPts, tchar, temp, 1, k[PBe9star], 1);//>
        nVar++;
        // Galerkin project solution back to original space
        // m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, temp, pb[PBe9star]);
        // for(int i = 0; i<nScaledDomainPts ; i++) 
        // {
        //     if(pb[PBe9star][i]<0) pb[PBe9star][i] = 0;
        // }

        //PBz9
        //N9b*Pmaxb*Pba/s9 + N9se*Pmaxse*Psea/s9
        Vmath::Svtsvtp(nScaledDomainPts, N9b*Pmaxb/s9, PLPba, 1, N9se*Pmaxse/s9, PLPsea, 1, temp, 1);
        //-ze9btot
        Vmath::Vsub(nScaledDomainPts, temp, 1, ze9btotInterp, 1, temp, 1);  
        //k9on*PBz9*v9*temp
        Vmath::Vmul(nScaledDomainPts, pb[PBz9], 1, temp, 1, temp2, 1);
        Vmath::Smul(nScaledDomainPts, -k9on*v9, temp2, 1, temp, 1);
        //-k9off*PBz9
        Vmath::Svtvp(nScaledDomainPts, -k9off, pb[PBz9], 1, temp, 1, k[PBz9], 1);
//        Vmath::Smul(nScaledDomainPts, tchar, temp, 1, k[PBz9], 1);//>
        nVar++;
        // Galerkin project solution back to original space
        // m_fields[0]->PhysGalerkinProjection1DScaled(OneDptscale, temp, pb[PBz9]);
        // for(int i = 0; i<nScaledDomainPts ; i++) 
        // {
        //     if(pb[PBz9][i]<0) pb[PBz9][i] = 0;
        // }

        //Subendothelium bound species:
        //Need bot row concentrations of certain fp and se species
        Array<OneD, Array<OneD, NekDouble> > BFP = Array<OneD, Array<OneD, NekDouble> > (nVariables);
        Array<OneD, Array<OneD, NekDouble> > BSE = Array<OneD, Array<OneD, NekDouble> > (numSEVars);

        for (int i = 0; i<nVariables; i++)
        {
            BFP[i] = Array<OneD, NekDouble>(nScaledBtr, 0.0);
            Vmath::Gathr(nScaledBtr, fpInterp[i], brCoord,  BFP[i]);
            // if(Vmath::Vmin(nScaledDomainPts, fp[i], 1) < 0) 
            //     cout<<"\n fp["<<i<<"<<] is negative inside BFP generation!";
            // if(Vmath::Vmin(nBtr, BFP[i],1)<0)
            //     cout<<"\n\n BFP["<<i<<"] is negative after BFP generation";
        }

        for (int i = 0; i<numSEVars; i++)
        {
            BSE[i] = Array<OneD, NekDouble>(nScaledBtr, 0.0);
            Vmath::Gathr(nScaledBtr, se[i], brCoord, BSE[i]);
        }
        
        
        //Note: SE chemicals only defined on bot row
        Array<OneD, NekDouble> tmp = Array<OneD, NekDouble>(nScaledBtr);
        Array<OneD, NekDouble> tmp2 = Array<OneD, NekDouble>(nScaledBtr);
        
        //SEz10e7 
        //kz10e7seplu FPz10 SEe7 - (kz10e7semin +kz10e7secat)SEz10e7
        Vmath::Vmul(nScaledBtr, BFP[FPz10], 1, BSE[SEe7], 1, tmp2, 1);
        Vmath::Svtsvtp(nScaledBtr, kz10e7seplu, tmp2, 1, -(kz10e7semin +kz10e7secat), BSE[SEz10e7], 1, tmp, 1);    
        
        //2D: Loop over quadrature points to assign max cover like:
        //cover(i) = kadh(i,1)*P0*max(Pma(i,1)+Pmu(i,1)+Pba(i,1)*(Pmaxb/P0), 
        //                            Pma(i,2)+Pmu(i,2)+Pba(i,2)*(Pmaxb/P0))
        
        for (int i = 0; i < nScaledBtr; i++)
        {   
            for (int j = 0; j< nScaledDomainPts; j++)
            {
                if(interpCoord_x[brCoord[i]] == interpCoord_x[j] 
                   && 
                   interpCoord_y[j] >= 0 
                   && 
                   interpCoord_y[j] <= m_rzt) 
                {
                    
                    if(cover[i] < fpInterp[Pma][j]+fpInterp[Pmu][j]+PLPba[j]*(Pmaxb/P0)) 
                        cover[i] = fpInterp[Pma][j]+fpInterp[Pmu][j]+PLPba[j]*(Pmaxb/P0);             
                }
            }
        }
        Vmath::Vmul(nScaledBtr, Bkadh, 1, cover, 1, cover, 1);
        Vmath::Smul(nScaledBtr, P0, cover, 1, cover, 1);
        
        //-kadh(x)SEz10e7 cover 
        Vmath::Vmul(nScaledBtr, Bkadh, 1, BSE[SEz10e7], 1, tmp2, 1);
        Vmath::Vmul(nScaledBtr, tmp2, 1, cover, 1, tmp2, 1);
        Vmath::Vsub(nScaledBtr, tmp, 1, tmp2, 1, tmp, 1);
        Vmath::Scatr(nScaledBtr, tmp, brCoord,  k[nVar]);
        nVar++;
//        Vmath::Smul(nScaledDomainPts, tchar, temp, 1, k[nVar++], 1);//>
        
        //SEe7 
        //tfavail = tf - z7se - e7se - z7se_e2 - z7se_e10  
        //          - z9_e7se - z10_e7se - tfpi_e10_e7se - e10_e7se
        Vmath::Vsub(nScaledBtr, BSE[SEtf], 1, BSE[SEz7], 1, tfAvail, 1);
        Vmath::Vsub(nScaledBtr, tfAvail, 1, BSE[SEe7], 1, tfAvail, 1);
        Vmath::Vsub(nScaledBtr, tfAvail, 1, BSE[SEz7e2], 1, tfAvail, 1);
        Vmath::Vsub(nScaledBtr, tfAvail, 1, BSE[SEz7e10], 1, tfAvail, 1);
        Vmath::Vsub(nScaledBtr, tfAvail, 1, BSE[SEz10e7], 1, tfAvail, 1);
        Vmath::Vsub(nScaledBtr, tfAvail, 1, BSE[SEtfpie10e7], 1, tfAvail, 1);
        //k7on*BFP[FPe7]*v7*tfavail  
        Vmath::Vmul(nScaledBtr, tfAvail, 1, BFP[FPe7], 1, tmp, 1);
        Vmath::Smul(nScaledBtr, k7on*v7, tmp, 1, tmp, 1);
        //-k7off*e7se + k2cat*z7se_e10
        Vmath::Svtsvtp(nScaledBtr, kz7e10secat, BSE[SEz7e10], 1, -k7off, 
                       BSE[SEe7], 1, tmp2, 1);
        Vmath::Vadd(nScaledBtr, tmp2, 1, tmp, 1, tmp, 1); 
        //+k3cat*z7se_e2
        Vmath::Svtvp(nScaledBtr, kz7e2secat, BSE[SEz7e2], 1, tmp, 1, tmp, 1);
        //+(k8min+k8cat)*z10_e7se
        Vmath::Svtvp(nScaledBtr, (kz10e7semin+kz10e7secat), BSE[SEz10e7], 
                     1, tmp, 1, tmp, 1);
        //-k8plu*BFP[FPz10]*v10*e7se 
        Vmath::Vmul(nScaledBtr, BFP[FPz10], 1, BSE[SEe7], 1, tmp2, 1);
        Vmath::Svtvp(nScaledBtr, -kz10e7seplu*v10, tmp2, 1, tmp, 1, tmp, 1);  
        //-k9plu*BFP[FPz9]*v9*e7se
        Vmath::Vmul(nScaledBtr, BFP[FPz9], 1, BSE[SEe7], 1, tmp2, 1);
        Vmath::Svtvp(nScaledBtr, -kz9e7seplu*v9, tmp2, 1, tmp, 1, tmp, 1) ; 
        //-k11plu*BFP[FPtfpie10]*vf*e7se
        Vmath::Vmul(nScaledBtr, BFP[FPtfpie10], 1, BSE[SEe7], 1, tmp2, 1);
        Vmath::Svtvp(nScaledBtr, -ktfpie10e7seplu*vf, tmp2, 1, tmp, 1, tmp, 1) ;         
        //+k11min*tfpi_e10_e7se
        Vmath::Svtvp(nScaledBtr, ktfpie10e7semin, BSE[SEtfpie10e7], 1, tmp, 1, tmp, 1) ;  
        //-e7se*cover
        Vmath::Vmul(nScaledBtr, BSE[SEe7], 1, cover, 1, tmp2, 1);
        Vmath::Vsub(nScaledBtr, tmp, 1, tmp2, 1, tmp, 1);

        Vmath::Scatr(nScaledBtr, tmp, brCoord,  k[nVar]);
        nVar++;

        //SEtfpie10e7
        //k11plu*BFP[FPtfpie10]*vf*e7se
        Vmath::Vmul(nScaledBtr, BFP[FPtfpie10], 1, BSE[SEe7], 1, tmp, 1);
        Vmath::Smul(nScaledBtr, ktfpie10e7seplu*vf, tmp, 1, tmp, 1);
        //-k11min*tfpi_e10_e7se 
        Vmath::Svtvp(nScaledBtr, -ktfpie10e7semin, BSE[SEtfpie10e7], 1, tmp, 1, tmp, 1) ;          
        //- tfpi_e10_e7se*cover
        Vmath::Vmul(nScaledBtr, BSE[SEtfpie10e7], 1, cover, 1, tmp2, 1);
        Vmath::Vsub(nScaledBtr, tmp, 1, tmp2, 1, tmp, 1);
        Vmath::Scatr(nScaledBtr, tmp, brCoord,  k[nVar]);
        nVar++;

        //SEz7
        //k7on*BFP[FPz7]*v7*tfavail
        Vmath::Vmul(nScaledBtr, tfAvail, 1, BFP[FPz7], 1, tmp, 1);

        //cout<<"\n max(tfAvail) = "<<Vmath::Vmax(nScaledBtr, tfAvail, 1)<< 
        //" max(BFP[FPz7]) = "<<Vmath::Vmax(nScaledBtr, BFP[FPz7], 1);
        Vmath::Smul(nScaledBtr, v7*k7on, tmp, 1, tmp, 1);
        //- k7off*z7se 
        Vmath::Svtvp(nScaledBtr, -k7off, BSE[SEz7], 1, tmp, 1, tmp, 1) ;          
        //- k2plu*z7se*BFP[FPe10]*v10 
        Vmath::Vmul(nScaledBtr, BSE[SEz7], 1, BFP[FPe10], 1, tmp2, 1);
        Vmath::Svtvp(nScaledBtr, -kz7e10seplu*v10, tmp2, 1, tmp, 1, tmp, 1) ;          
        //+k2min*z7se_e10 
        Vmath::Svtvp(nScaledBtr, kz7e10semin, BSE[SEz7e10], 1, tmp, 1, tmp, 1) ;          
        //- k3plu*z7se*BFP[FPe2]*v2 
        Vmath::Vmul(nScaledBtr, BSE[SEz7], 1, BFP[FPe2], 1, tmp2, 1);
        Vmath::Svtvp(nScaledBtr, -kz7e2seplu*v2, tmp2, 1, tmp, 1, tmp, 1) ;          
        //+ k3min*z7se_e2 
        Vmath::Svtvp(nScaledBtr, kz7e2semin, BSE[SEz7e2], 1, tmp, 1, tmp, 1) ;          
        //- z7se*cover
        Vmath::Vmul(nScaledBtr, BSE[SEz7], 1, cover, 1, tmp2, 1);
        Vmath::Vsub(nScaledBtr, tmp, 1, tmp2, 1, tmp, 1);
        Vmath::Scatr(nScaledBtr, tmp, brCoord, k[nVar]);
        nVar++;

        //SEz7e10
        //k2plu*z7se*BFP[FPe10]*v10
        Vmath::Vmul(nScaledBtr, BFP[FPe10], 1, BSE[SEz7], 1, tmp, 1);
        Vmath::Smul(nScaledBtr, kz7e10seplu*v10, tmp, 1, tmp, 1);
        //-(k2min+k2cat)*z7se_e10 
        Vmath::Svtvp(nScaledBtr, -(kz7e10semin+kz7e10secat), 
                     BSE[SEz7e10], 1, tmp, 1, tmp, 1) ;          
        //- z7se_e10*cover 
        Vmath::Vmul(nScaledBtr, cover, 1, BSE[SEz7e10], 1, tmp2, 1);
        Vmath::Vsub(nScaledBtr, tmp, 1, tmp2, 1, tmp, 1);
        Vmath::Scatr(nScaledBtr, tmp, brCoord,  k[nVar]);
        nVar++;

        //SEtf
        //-tf*cover
        Vmath::Vmul(nScaledBtr, BSE[SEtf], 1, cover, 1, tmp, 1);
        Vmath::Neg(nScaledBtr, tmp, 1);
        Vmath::Scatr(nScaledBtr, tmp, brCoord,  k[nVar]);
        nVar++;

        //SEz7e2
        //k3plu*z7se*BFP[FPe2]*v2
        Vmath::Vmul(nScaledBtr, BSE[SEz7], 1, BFP[FPe2], 1, tmp, 1);
        Vmath::Smul(nScaledBtr, kz7e2seplu*v2, tmp, 1, tmp, 1);
        //- (k3min+k3cat)*z7se_e2 
        Vmath::Svtvp(nScaledBtr, -(kz7e2semin+kz7e2secat), BSE[SEz7e2], 
                     1, tmp, 1, tmp, 1);
        //- z7se_e2*cover
        Vmath::Vmul(nScaledBtr, BSE[SEz7e2], 1, cover, 1, tmp2, 1);
        Vmath::Vsub(nScaledBtr, tmp, 1, tmp2, 1, tmp, 1);
        Vmath::Scatr(nScaledBtr, tmp, brCoord,  k[nVar]);
        nVar++;
        
        //SEz9e7
        //k9plu*e7se*Bz9*v9 
        Vmath::Vmul(nScaledBtr, BSE[SEe7], 1, BFP[FPz9], 1, tmp, 1);
        Vmath::Smul(nScaledBtr, kz9e7seplu*v9, tmp, 1, tmp, 1);
        //-(k9min+k9cat)*z9_e7se 
        Vmath::Svtvp(nScaledBtr, -(kz9e7semin+kz9e7secat), BSE[SEz9e7], 
                     1, tmp, 1, tmp, 1);
        //-z9_e7se*cover
        Vmath::Vmul(nScaledBtr, cover, 1, BSE[SEz9e7], 1, tmp2, 1);
        Vmath::Vsub(nScaledBtr, tmp, 1, tmp2, 1, tmp, 1);
        Vmath::Scatr(nScaledBtr, tmp, brCoord,  k[nVar]);
        nVar++;
        
        // 4. Pba = kcoh*P0/Pmaxse*geta*Pma - kadh*Pmaxse*Pmat2*Pba
        Vmath::Vmul(nScaledDomainPts, fpInterp[Pma], 1, geta, 1, temp2, 1); // Pma*geta
        Vmath::Smul(nScaledDomainPts, kcoh*P0/Pmaxse, temp2, 1, temp2, 1);
        Vmath::Sadd(nScaledDomainPts, -1.0, PLPsea, 1, temp, 1); //Psea-1
        Vmath::Vmul(nScaledDomainPts, temp, 1, PLPba, 1, temp, 1); //Pba(Psea-1)
        Vmath::Vmul(nScaledDomainPts, temp, 1, kadh, 1, temp, 1); //kadh Pba (Psea-1)
        Vmath::Svtvp(nScaledDomainPts, Pmaxse, temp, 1, temp2, 1,  k[nVar], 1);
        nVar++;
	
        // 3. Psea = kadh*P0*(1-Psea)*(Pma+Pmu+Pba*(Pmaxse/P0))
        Vmath::Vadd(nScaledDomainPts, fpInterp[Pma], 1, fpInterp[Pmu], 1 , temp, 1); //Pmu+Pma
        Vmath::Smul(nScaledDomainPts, Pmaxse/P0, PLPba, 1, temp2, 1);
        Vmath::Vadd(nScaledDomainPts, temp, 1, temp2, 1, temp, 1); // temp = Pmu+Pma+Pba*(Pmax/P0)
        
        Vmath::Sadd(nScaledDomainPts, -1.0, PLPsea, 1, temp2, 1); // Psea-1
        Vmath::Vmul(nScaledDomainPts, kadh, 1, temp2, 1, temp2, 1); //temp2 = kadh(Psea-1)
        Vmath::Smul(nScaledDomainPts, P0, temp2, 1, temp2, 1); //kadh(Psea-1)P0
        Vmath::Neg(nScaledDomainPts, temp2, 1); //kadh(1-Psea)P0
        Vmath::Vmul(nScaledDomainPts, temp2, 1, temp, 1,  k[nVar], 1); //done
        
        //check if any k is less that 1e-14, then set it to 0
        for(int i = 0; i <numSEVars+numPBVars+2; i++)
            for(int j = 0; j<nScaledDomainPts; j++)
                if(abs(k[i][j])<1e-14) 
                    k[i][j] = 0;
    }
  
/**
 * @Brief Update Release term in reaction part of ADP
 * 
 * @param inarray    Given fields.
 * @param sigmaz     Output: current release value
 * @param sigmar     Holds future info
 */ 
    void UnsteadyAdvectionDiffusion::UpdateADPRelease()
    {
        /*            FOR ADP RELEASE
         *     length of Rtemp  -    5/rlen = 0.25
         * this means each time to release is .25 seconds
         */
        //Array<OneD, NekDouble> PbnewInterp = Array<OneD, NekDouble>(nScaledDomainPts);
        //m_fields[0]->PhysInterp1DScaled(OneDptscale, Pbnew, PbnewInterp);

        for( int i = 0; i < nScaledDomainPts; i++)
        {
            Array<OneD, NekDouble> sigmatemp(rlen);
            // sigmatemp(1:rlen-1)=sigmar(i,j,2:rlen)
            Vmath::Vcopy(rlen-1, &sigmar[i][1], 1, &sigmatemp[0], 1);
            // sigmar(i,j,:)= Pbnew(i,j)*rtemp*Pmaxse*ADPrel*6.022e23 + sigmatemp
            Vmath::Svtvp(rlen, Pbnew[i]*Pmaxse*ADPrel*6.022e23, rtemp, 1, sigmatemp, 1, sigmar[i], 1);
            // sigmaz(i,j)=sigmar(i,j,1)
            //	     Vmath::Vcopy(rlen, &sigmar[i][0], 1, &sigmaz[i], 1);
            sigmaz[i] = sigmar[i][0];
        }
        int nSolutionPts = GetNpoints();
        Pbnew = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
    }

/**
 * @brief Compute the projection for the unsteady advection 
 * diffusion problem.
 * 
 * @param inarray    Given fields.
 * @param outarray   Calculated solution.
 * @param time       Time.
 */
    void UnsteadyAdvectionDiffusion::DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble time)
    {
        int i;
        int nvariables = inarray.num_elements();
        PACSetBoundaryConditions();
        //SetBoundaryConditions(time);
        switch(m_projectionType)
        {
        case MultiRegions::eDiscontinuous:
        {
            // Just copy over array
            int npoints = GetNpoints();

            for(i = 0; i < nvariables; ++i)
            {
                Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
            }
            break;
        }
        case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
        {
            Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());

            for(i = 0; i < nvariables; ++i)
            {
                m_fields[i]->FwdTrans(inarray[i], coeffs);
                m_fields[i]->BwdTrans_IterPerExp(coeffs, outarray[i]);
            }
            break;
        }
        default:
        {
            ASSERTL0(false, "Unknown projection scheme");
            break;
        }
        }
    }

    /**
     * If boundary conditions are field-dependent, they will be evaluated here
     */
    void UnsteadyAdvectionDiffusion::PACSetBoundaryConditions()
    {
        //std::string varName;
        int nvariables = m_fields.num_elements();
        for (int i = 0; i < nvariables; ++i)
        {
            Array<OneD, NekDouble> cVal, cValOutFlow; 
            FillRobinBoundaryCoefficients(i, cVal, cValOutFlow);
            m_fields[i]->PACEvaluateBoundaryConditions(i, cVal, cValOutFlow);
            
        }
    }

    /**
     * Robin boundary coefficients are values of certain SE and FP chemicals
     * depending on which variable we are applying robin boundary to.  
     * @param varID  variable ID whose robin boundary condition is required
     */
    void UnsteadyAdvectionDiffusion::FillRobinBoundaryCoefficients(int varID, Array<OneD, NekDouble> &cVal, Array<OneD, NekDouble> &cValOutFlow)
    {
        //cout<<"\n Inside fillrobinboundaryconditions: m_epsilon = "<<m_epsilon;
        switch(varID)
        {
        case Pmu:
        case Pma:{
            cValOutFlow = Array<OneD, NekDouble>(nScOutflow);
            Array<OneD, NekDouble> temp = Array<OneD, NekDouble>(nScOutflow);
            Array<OneD, NekDouble> mobileAdvVel_x = Array<OneD, NekDouble>(nSolutionPts, nonDAdvVel);
            Array<OneD, NekDouble> mobileAdvVelInt_x = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
            // velo = W(\phi_T)u
            Vmath::Vmul(nSolutionPts, m_velocity[0], 1, mobileAdvVel_x, 1, mobileAdvVel_x, 1);
            Vmath::Vmul(nSolutionPts, phiT, 1, mobileAdvVel_x, 1, mobileAdvVel_x, 1);
            m_fields[0]->PhysInterp1DScaled(OneDptscale, mobileAdvVel_x, mobileAdvVelInt_x);

            Vmath::Gathr(nScOutflow, mobileAdvVelInt_x, ofCoord, temp);
            Vmath::Gathr(nScOutflow, fpInterp[Pmu], ofCoord, cValOutFlow);
            Vmath::Vmul(nScOutflow, temp, 1, cValOutFlow, 1, cValOutFlow, 1);
            //Vmath::Neg(nScOutflow, cValOutFlow, 1);
            break;}
        case FPe2:  {// ref eq A.51 left to right
            cVal = Array<OneD, NekDouble>(numRZ);
            Array<OneD, NekDouble> temp = Array<OneD, NekDouble>(numRZ);
            Vmath::Gathr(numRZ, se[SEz7], rzCoord, temp); //z7se
             //cout<<"\n Robin boundary SEz7 max = "<<Vmath::Vmax(numRZ, temp, 1);
            Vmath::Gathr(numRZ, fpInterp[FPe2], rzCoord, cVal); //e2
            //cout<<"\t Robin boundary FPe2 max = "<<Vmath::Vmax(numRZ, cVal, 1);
            Vmath::Vmul(numRZ, temp, 1, cVal, 1, cVal, 1);
            Vmath::Smul(numRZ, s7d*kz7e2seplu/m_epsilon, cVal, 1, cVal, 1);
            Vmath::Gathr(numRZ, se[SEz7e2], rzCoord, temp); //z7e2se
            Vmath::Svtvp(numRZ, -(s7d/v2)*(kz7e2semin + kz7e2secat)/m_epsilon, se[SEz7e2], 1, cVal, 1, cVal, 1);
            //if(Vmath::Vmin(numRZ, cVal, 1) < 0)
             //   cout<<"\t Final  cVal max = "<<Vmath::Vmax(numRZ, cVal, 1);
            //Vmath::Neg(numRZ, cVal, 1);
            break;}
        case FPz5:{
            cValOutFlow = Array<OneD, NekDouble>(nScOutflow);
            Array<OneD, NekDouble> temp = Array<OneD, NekDouble>(nScOutflow);
            Array<OneD, NekDouble> mobileAdvVel_x = Array<OneD, NekDouble>(nSolutionPts, nonDAdvVel);
            Array<OneD, NekDouble> mobileAdvVelInt_x = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
            Vmath::Vmul(nSolutionPts, m_velocity[0], 1, mobileAdvVel_x, 1, mobileAdvVel_x, 1);
            m_fields[0]->PhysInterp1DScaled(OneDptscale, mobileAdvVel_x, mobileAdvVelInt_x);
            Vmath::Gathr(nScOutflow, mobileAdvVelInt_x, ofCoord, temp);
            Vmath::Gathr(nScOutflow, fpInterp[FPz5], ofCoord, cValOutFlow);
            Vmath::Vmul(nScOutflow, temp, 1, cValOutFlow, 1, cValOutFlow, 1);
           // Vmath::Neg(nScOutflow, cValOutFlow, 1);
            break;}
        case FPtfpi:{
            cValOutFlow = Array<OneD, NekDouble>(nScOutflow);
            Array<OneD, NekDouble> temp = Array<OneD, NekDouble>(nScOutflow);
            Array<OneD, NekDouble> mobileAdvVel_x = Array<OneD, NekDouble>(nSolutionPts, nonDAdvVel);
            Array<OneD, NekDouble> mobileAdvVelInt_x = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
            Vmath::Vmul(nSolutionPts, m_velocity[0], 1, mobileAdvVel_x, 1, mobileAdvVel_x, 1);
            m_fields[0]->PhysInterp1DScaled(OneDptscale, mobileAdvVel_x, mobileAdvVelInt_x);
            Vmath::Gathr(nScOutflow, mobileAdvVelInt_x, ofCoord, temp);
            Vmath::Gathr(nScOutflow, fpInterp[FPtfpi], ofCoord, cValOutFlow);
            Vmath::Vmul(nScOutflow, temp, 1, cValOutFlow, 1, cValOutFlow, 1);
            //Vmath::Neg(nScOutflow, cValOutFlow, 1);
            break;}
        case FPz2:{
            cValOutFlow = Array<OneD, NekDouble>(nScOutflow);
            Array<OneD, NekDouble> temp = Array<OneD, NekDouble>(nScOutflow);
            Array<OneD, NekDouble> mobileAdvVel_x = Array<OneD, NekDouble>(nSolutionPts, nonDAdvVel);
            Array<OneD, NekDouble> mobileAdvVelInt_x = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
            Vmath::Vmul(nSolutionPts, m_velocity[0], 1, mobileAdvVel_x, 1, mobileAdvVel_x, 1);
            m_fields[0]->PhysInterp1DScaled(OneDptscale, mobileAdvVel_x, mobileAdvVelInt_x);
            Vmath::Gathr(nScOutflow, mobileAdvVelInt_x, ofCoord, temp);
            Vmath::Gathr(nScOutflow, fpInterp[FPz2], ofCoord, cValOutFlow);
            Vmath::Vmul(nScOutflow, temp, 1, cValOutFlow, 1, cValOutFlow, 1);
           // Vmath::Neg(nScOutflow, cValOutFlow, 1);
            break;}
        case FPz7:  {// ref eq A.52 left to right
            cVal = Array<OneD, NekDouble>(numRZ);
            Array<OneD, NekDouble> temp = Array<OneD, NekDouble>(numRZ);
            Array<OneD, NekDouble> scaledA = Array<OneD, NekDouble>(nScaledDomainPts);
            Vmath::Gathr(numRZ, fpInterp[FPz7], rzCoord, temp); //z7

            //first, scatter tfAvail on large grid then extract only RZ
            //Vmath::Scatr(nBtr, tfAvail, brCoord, scaledA);
            Vmath::Gathr(numRZ, tfAvail, rzCoord, cVal); //tfAvail
                        
            Vmath::Vmul(numRZ, temp, 1, cVal, 1, cVal, 1);
            Vmath::Smul(numRZ, s7d*k7on/m_epsilon, cVal, 1, cVal, 1);
            Vmath::Gathr(numRZ, se[SEz7], rzCoord, temp); //z7se
            Vmath::Svtvp(numRZ, -(k7off/m_epsilon)*(s7d/v7), temp, 1, cVal, 1, cVal, 1);
            //Vmath::Neg(numRZ, cVal, 1);
           
            
            cValOutFlow = Array<OneD, NekDouble>(nScOutflow);
            Array<OneD, NekDouble> temp2 = Array<OneD, NekDouble>(nScOutflow);
            Array<OneD, NekDouble> mobileAdvVel_x = Array<OneD, NekDouble>(nSolutionPts, nonDAdvVel);
            Array<OneD, NekDouble> mobileAdvVelInt_x = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
            Vmath::Vmul(nSolutionPts, m_velocity[0], 1, mobileAdvVel_x, 1, mobileAdvVel_x, 1);
            m_fields[0]->PhysInterp1DScaled(OneDptscale, mobileAdvVel_x, mobileAdvVelInt_x);
            Vmath::Gathr(nScOutflow, mobileAdvVelInt_x, ofCoord, temp2);
            Vmath::Gathr(nScOutflow, fpInterp[FPz7], ofCoord, cValOutFlow);
            Vmath::Vmul(nScOutflow, temp2, 1, cValOutFlow, 1, cValOutFlow, 1);
            //Vmath::Neg(nScOutflow, cValOutFlow, 1);
            break;}
        case FPz8:{
            cValOutFlow = Array<OneD, NekDouble>(nScOutflow);
            Array<OneD, NekDouble> temp = Array<OneD, NekDouble>(nScOutflow);
            Array<OneD, NekDouble> mobileAdvVel_x = Array<OneD, NekDouble>(nSolutionPts, nonDAdvVel);
            Array<OneD, NekDouble> mobileAdvVelInt_x = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
            Vmath::Vmul(nSolutionPts, m_velocity[0], 1, mobileAdvVel_x, 1, mobileAdvVel_x, 1);
            m_fields[0]->PhysInterp1DScaled(OneDptscale, mobileAdvVel_x, mobileAdvVelInt_x);
            Vmath::Gathr(nScOutflow, mobileAdvVelInt_x, ofCoord, temp);
            Vmath::Gathr(nScOutflow, fpInterp[FPz8], ofCoord, cValOutFlow);
            Vmath::Vmul(nScOutflow, temp, 1, cValOutFlow, 1, cValOutFlow, 1);
           // Vmath::Neg(nScOutflow, cValOutFlow, 1);
            break;}
        case FPe7:
        {// ref eq A.53 left to right
            cVal = Array<OneD, NekDouble>(numRZ);
            Array<OneD, NekDouble> temp = Array<OneD, NekDouble>(numRZ);
            Array<OneD, NekDouble> scaledA = Array<OneD, NekDouble>(nScaledDomainPts);
            Vmath::Gathr(numRZ, fpInterp[FPe7], rzCoord, temp); //e7

            //first, scatter tfAvail on large grid then extract only RZ
            //Vmath::Scatr(nScaledBtr, tfAvail, brCoord, scaledA);
            // no need to scatr coz tfAvail is already defined on brCoord
            Vmath::Gathr(numRZ, tfAvail, rzCoord, cVal); //tfAvail
                        
            Vmath::Vmul(numRZ, temp, 1, cVal, 1, cVal, 1);
            Vmath::Smul(numRZ, s7d*k7on/m_epsilon, cVal, 1, cVal, 1);
            Vmath::Gathr(numRZ, se[SEz7], rzCoord, temp); //z7se
            Vmath::Svtvp(numRZ, -(k7off/m_epsilon)*(s7d/v7), temp, 1, cVal, 1, cVal, 1);
            //Vmath::Neg(numRZ, cVal, 1);
           

           cValOutFlow = Array<OneD, NekDouble>(nScOutflow);
            Array<OneD, NekDouble> temp2 = Array<OneD, NekDouble>(nScOutflow);
            Array<OneD, NekDouble> mobileAdvVel_x = Array<OneD, NekDouble>(nSolutionPts, nonDAdvVel);
            Array<OneD, NekDouble> mobileAdvVelInt_x = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
            Vmath::Vmul(nSolutionPts, m_velocity[0], 1, mobileAdvVel_x, 1, mobileAdvVel_x, 1);
            m_fields[0]->PhysInterp1DScaled(OneDptscale, mobileAdvVel_x, mobileAdvVelInt_x);
            Vmath::Gathr(nScOutflow, mobileAdvVelInt_x, ofCoord, temp2);
            Vmath::Gathr(nScOutflow, fpInterp[FPe7], ofCoord, cValOutFlow);
            Vmath::Vmul(nScOutflow, temp2, 1, cValOutFlow, 1, cValOutFlow, 1);
           // Vmath::Neg(nScOutflow, cValOutFlow, 1);
            break;}

        case FPz9:  {// ref eq A.54 left to right
            cVal = Array<OneD, NekDouble>(numRZ);
            Array<OneD, NekDouble> temp = Array<OneD, NekDouble>(numRZ);
            Vmath::Gathr(numRZ, fpInterp[FPz7], rzCoord, cVal); //z9
            Vmath::Gathr(numRZ, se[SEe7], rzCoord, temp); //z7se
            Vmath::Vmul(numRZ, temp, 1, cVal, 1, cVal, 1);
            Vmath::Smul(numRZ, s7d*kz9e7seplu/m_epsilon, cVal, 1, cVal, 1);

            Vmath::Gathr(numRZ, se[SEz9e7], rzCoord, temp); //z9e7se
            Vmath::Svtvp(numRZ, -(kz9e7semin/m_epsilon)*(s7d/v9), temp, 1, cVal, 1, cVal, 1);
           // Vmath::Neg(numRZ, cVal, 1);
           
 
            cValOutFlow = Array<OneD, NekDouble>(nScOutflow);
            Array<OneD, NekDouble> temp2 = Array<OneD, NekDouble>(nScOutflow);
            Array<OneD, NekDouble> mobileAdvVel_x = Array<OneD, NekDouble>(nSolutionPts, nonDAdvVel);
            Array<OneD, NekDouble> mobileAdvVelInt_x = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
            Vmath::Vmul(nSolutionPts, m_velocity[0], 1, mobileAdvVel_x, 1, mobileAdvVel_x, 1);
            m_fields[0]->PhysInterp1DScaled(OneDptscale, mobileAdvVel_x, mobileAdvVelInt_x);
            Vmath::Gathr(nScOutflow, mobileAdvVelInt_x, ofCoord, temp2);
            Vmath::Gathr(nScOutflow, fpInterp[FPz9], ofCoord, cValOutFlow);
            Vmath::Vmul(nScOutflow, temp2, 1, cValOutFlow, 1, cValOutFlow, 1);
           // Vmath::Neg(nScOutflow, cValOutFlow, 1);
            break;}
        case FPe9:   {// ref eq A.55 left to right
            cVal = Array<OneD, NekDouble>(numRZ);
            Vmath::Gathr(numRZ, se[SEz9e7], rzCoord, cVal); //z9e7se
            Vmath::Smul(numRZ, -(kz9e7secat/m_epsilon)*s7d/v9, cVal, 1, cVal, 1);
            break;}
        case FPz10:  {// ref eq A.56 left to right
            cVal = Array<OneD, NekDouble>(numRZ);
            Array<OneD, NekDouble> temp = Array<OneD, NekDouble>(numRZ);
            Vmath::Gathr(numRZ, se[SEe7], rzCoord, temp); //e7se
            Vmath::Gathr(numRZ, fpInterp[FPz10], rzCoord, cVal); //z10
            Vmath::Vmul(numRZ, temp, 1, cVal, 1, cVal, 1);
            Vmath::Smul(numRZ, s7d*(kz10e7seplu/m_epsilon), cVal, 1, cVal, 1);

            Vmath::Gathr(numRZ, se[SEz10e7], rzCoord, temp); //z10e7se
            Vmath::Svtvp(numRZ, -(s7d/v10)*(kz10e7semin/m_epsilon), temp, 1, cVal, 1, cVal, 1);
            //Vmath::Neg(numRZ, cVal, 1);
           

            cValOutFlow = Array<OneD, NekDouble>(nScOutflow);
            Array<OneD, NekDouble> temp2 = Array<OneD, NekDouble>(nScOutflow);
            Array<OneD, NekDouble> mobileAdvVel_x = Array<OneD, NekDouble>(nSolutionPts, nonDAdvVel);
            Array<OneD, NekDouble> mobileAdvVelInt_x = Array<OneD, NekDouble>(nScaledDomainPts, 0.0);
            Vmath::Vmul(nSolutionPts, m_velocity[0], 1, mobileAdvVel_x, 1, mobileAdvVel_x, 1);
            m_fields[0]->PhysInterp1DScaled(OneDptscale, mobileAdvVel_x, mobileAdvVelInt_x);
            Vmath::Gathr(nScOutflow, mobileAdvVelInt_x, ofCoord, temp2);
            Vmath::Gathr(nScOutflow, fpInterp[FPz10], ofCoord, cValOutFlow);
            Vmath::Vmul(nScOutflow, temp2, 1, cValOutFlow, 1, cValOutFlow, 1);
           // Vmath::Neg(nScOutflow, cValOutFlow, 1);
           
            break;}
        case FPe10:  {// ref eq A.57 left to right
            cVal = Array<OneD, NekDouble>(numRZ);
            Array<OneD, NekDouble> temp = Array<OneD, NekDouble>(numRZ);
            Vmath::Gathr(numRZ, se[SEe7], rzCoord, temp); //SEz7
            Vmath::Gathr(numRZ, fpInterp[FPe10], rzCoord, cVal); //e10
            Vmath::Vmul(numRZ, temp, 1, cVal, 1, cVal, 1);
            Vmath::Smul(numRZ, s7d*(kz7e10seplu/m_epsilon), cVal, 1, cVal, 1);
            
            Vmath::Gathr(numRZ, se[SEz10e7], rzCoord, temp); //z10e7se
            Vmath::Svtvp(numRZ, -(s7d/v10)*(kz10e7secat/m_epsilon), temp, 1, cVal, 1, cVal, 1);

            Vmath::Gathr(numRZ, se[SEz7e10], rzCoord, temp); //z7e10se
            Vmath::Svtvp(numRZ, -(s7d/v10)*((kz7e10semin + kz7e10secat)/m_epsilon), temp, 1, cVal, 1, cVal, 1);
            //Vmath::Neg(numRZ, cVal, 1);
           
            break;}
        case FPtfpie10:
        {// ref eq A.58 left to right
            cVal = Array<OneD, NekDouble>(numRZ);
            Array<OneD, NekDouble> temp = Array<OneD, NekDouble>(numRZ);
            Vmath::Gathr(numRZ, se[SEe7], rzCoord, temp); //e7se
            Vmath::Gathr(numRZ, fpInterp[FPtfpie10], rzCoord, cVal); //FPtfpie10
            Vmath::Vmul(numRZ, temp, 1, cVal, 1, cVal, 1);
            Vmath::Smul(numRZ, s7d*(ktfpie10e7seplu/m_epsilon), cVal, 1, cVal, 1);

            Vmath::Gathr(numRZ, se[SEtfpie10e7], rzCoord, temp); //tfpie10e7se
            Vmath::Svtvp(numRZ, -(s7d/vf)*(ktfpie10e7semin/m_epsilon), temp, 1, cVal, 1, cVal, 1);
          //  Vmath::Neg(numRZ, cVal, 1);
           
            break;}
            
        default:    break;
        }
    }

    
/* @brief Compute the diffusion term implicitly. 
 * 
 * @param inarray    Given fields.
 * @param outarray   Calculated solution.
 * @param time       Time.
 * @param lambda     Diffusion coefficient.
 */
    void UnsteadyAdvectionDiffusion::DoImplicitSolve(
        const Array<OneD, const Array<OneD, NekDouble> >&inarray,
        Array<OneD,       Array<OneD, NekDouble> >&outarray,
        const NekDouble time,
        const NekDouble lambda)
    {

        //cout<<"\n******DOIMPLICITSOLVE*************";
        int nVariables = inarray.num_elements(); // only Pmu,Pma,Peta,ADP and FPe2
        double xchar, tchar, nonDDiffCoeff, nonDLambda;
        // Number of solution points
        int nSolutionPts = GetNpoints();

        // for(int i = 0; i<inarray.num_elements(); i++)
        // for(int k = 0; k<nSolutionPts ; k++) // fpReactInterp[Pmu] = max(0,fpReactInterp[Pmu])
        // {
        //     if(inarray[i][k]<0) {cout<<"\n Diffusion got negative val for i = "<<i;exit(0);};
        // }

	      
        Array<OneD, Array<OneD, NekDouble> > phiT(2); 
        Array<OneD, Array<OneD, NekDouble> > phiT1(2);
	
        StdRegions::ConstFactorMap factors;
	
        xchar = m_session->GetParameter("xchar");
        tchar = m_session->GetParameter("tchar");
        //cout<<"\n Inside DoImplicitSolve : m_epsilon = "<<m_epsilon; 
        nonDDiffCoeff = tchar*m_epsilon/(pow(xchar,2)); //ND-DiffusionCoefficient
        nonDLambda = lambda;//lambda/tchar;
        if(m_useSpecVanVisc)
        {
            factors[StdRegions::eFactorSVVCutoffRatio] = m_sVVCutoffRatio;
            factors[StdRegions::eFactorSVVDiffCoeff]   = m_sVVDiffCoeff/m_epsilon;
        }

        Array<OneD, Array< OneD, NekDouble> > F(nVariables);
        F[0] = Array<OneD, NekDouble> (nSolutionPts*nVariables);
        
        for (int n = 1; n < nVariables; ++n)
        {
            F[n] = F[n-1] + nSolutionPts;
        }
        
        //Setting boundary conditions
        PACSetBoundaryConditions();

        // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
        // inarray = input: \hat{rhs} -> output: \hat{Y}
        // outarray = output: nabla^2 \hat{Y}
        // where \hat = modal coeffs
        for (int i = 0; i < nVariables; ++i)
        {
            switch(i)
            {
            case Peta: 
            {    nonDDiffCoeff = tchar*m_epsilonPeta*m_epsilon/(pow(xchar,2)); //ND-DiffusionCoefficient
                    factors[StdRegions::eFactorLambda] = 1.0/nonDLambda/nonDDiffCoeff;
                    Vmath::Smul(nSolutionPts, -factors[StdRegions::eFactorLambda], 
                                inarray[i], 1, F[i], 1);
                    m_fields[Peta]->HelmSolve(F[Peta], m_fields[Peta]->UpdateCoeffs(), 
                                              NullFlagList, factors);
                    m_fields[Peta]->BwdTrans(m_fields[Peta]->GetCoeffs(), outarray[i]);
                    break;
            }
            case Pmu:
            case Pma:
            {    
                phiT[i] = Array<OneD, NekDouble> (nSolutionPts, 0.0);
                phiT1[i] = Array<OneD, NekDouble> (nSolutionPts, 0.0);
                Vmath::Vadd(nSolutionPts, inarray[Pmu], 1, 
                            inarray[Pma], 1, phiT[i], 1);
                Vmath::Smul(nSolutionPts, phitfac, 
                            phiT[i], 1, phiT[i], 1);
                Vmath::Vadd(nSolutionPts, PLPsea, 1, 
                            phiT[i], 1, phiT[i], 1);
                Vmath::Vadd(nSolutionPts, PLPba, 1, 
                            phiT[i], 1, phiT[i], 1);
                Vmath::Neg(nSolutionPts, phiT[i], 1);
                Vmath::Sadd(nSolutionPts, 1.0, phiT[i], 1, 
                            phiT[i], 1); /* 1-phiT*/
                Vmath::Smul(nSolutionPts, PI, 
                            phiT[i], 1, phiT[i], 1); /*PI*(1-phiT)*/
              
                transform(phiT[i].begin(), phiT[i].end(), phiT[i].begin(), NekTanh);/*WphiT*/
                factors[StdRegions::eFactorLambda] = 1.0/nonDLambda/nonDDiffCoeff;
        
                Vmath::Smul(nSolutionPts, 
                            -factors[StdRegions::eFactorLambda], 
                            inarray[i], 1, 
                            m_fields[i]->UpdatePhys(), 1);
                StdRegions::VarCoeffMap m_varcoeff;
                m_varcoeff[StdRegions::eVarCoeffD00] = phiT[i];
                m_varcoeff[StdRegions::eVarCoeffD11] = phiT[i];
                m_varcoeff[StdRegions::eVarCoeffD01] = Array<OneD, NekDouble>(nSolutionPts, 0.0);
                NekDouble cpuTime = 0.0, elapsed = 0.0; 
                Timer timer;
              
                //timer.Start();
                m_fields[i]->HelmSolve(m_fields[i]->GetPhys(),
                                       m_fields[i]->UpdateCoeffs(), 
                                       NullFlagList, 
                                       factors, 
                                       m_varcoeff);
                /*timer.Stop();
                  elapsed  = timer.TimePerTest(1);
                  cpuTime += elapsed;
                  cout << "-------------------------------------------" << endl;
                  cout << "Total Computation Time for Pmu or Pma  = " << cpuTime << "s" << endl;
                  cout << "-------------------------------------------" << endl;*/
               
                m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), 
                                      outarray[i]);//m_fields[i]->UpdatePhys());
                    
                //m_fields[i]->SetPhysState(false);
                    
                //outarray[i] = m_fields[i]->GetPhys();
                break;
            }
            case ADP:
            {
                nonDDiffCoeff = tchar*m_epsilonADP*m_epsilon/(pow(xchar,2)); //ND-DiffusionCoefficient 
                factors[StdRegions::eFactorLambda] = 1.0/nonDLambda/nonDDiffCoeff;
                Vmath::Smul(nSolutionPts, -factors[StdRegions::eFactorLambda], 
                            inarray[i], 1, F[i], 1);

                NekDouble cpuTime = 0.0, elapsed = 0.0; 
                Timer timer;

                //timer.Start();
                m_fields[ADP]->HelmSolve(F[ADP], m_fields[ADP]->UpdateCoeffs(), 
                                         NullFlagList, factors);

                /*       timer.Stop();

                         elapsed  = timer.TimePerTest(1);
                         cpuTime += elapsed;

                         cout << "-------------------------------------------" << endl;
                         cout << "Total Computation Time for ADP  = " << cpuTime << "s" << endl;
                         cout << "-------------------------------------------" << endl;*/
                m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), outarray[i]);
                break;
            }
            
            default: // all other FP chemicals
            {    nonDDiffCoeff = tchar*m_epsilon/(pow(xchar,2)); //ND-DiffusionCoefficient
                factors[StdRegions::eFactorLambda] = 1.0/nonDLambda/nonDDiffCoeff;
                Vmath::Smul(nSolutionPts, -factors[StdRegions::eFactorLambda], 
                            inarray[i], 1, F[i], 1);
                m_fields[i]->HelmSolve(F[i], m_fields[i]->UpdateCoeffs(), 
                                       NullFlagList, factors);
                m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), outarray[i]);
                break;

            }
            }

            // for(int k = 0; k<nSolutionPts ; k++) // fpReactInterp[Pmu] = max(0,fpReactInterp[Pmu])
            // {
            //     if(outarray[i][k]<0) {outarray[i][k] = 0;cout<<"\n Diffusion solve returned negative val for i = "<<i;};
            // }

        }
    
        //cout<<"\nDone Diffusion. Max and Mins are:";
        //cout<<"pb[PBz8]: "<< Vmath::Vmax(nSolutionPts, pb[PBz8],1)<<" Min ="<<Vmath::Vmin(nSolutionPts, pb[PBz8],1)<<"\t";
        //cout<<"FPe2: "<< Vmath::Vmax(nSolutionPts, outarray[FPe2],1)<<" Min ="<<Vmath::Vmin(nSolutionPts, outarray[FPe2],1)<<"\t";
        //cout<<"Pmu: "<< Vmath::Vmax(nSolutionPts, outarray[Pmu],1)<<" Min ="<<Vmath::Vmin(nSolutionPts, outarray[Pmu],1)<<"\t";
        //cout<<"Pma: "<< Vmath::Vmax(nSolutionPts, outarray[Pma],1)<<" Min ="<<Vmath::Vmin(nSolutionPts, outarray[Pma],1)<<"\t";
        //cout<<"Peta: "<< Vmath::Vmax(nSolutionPts, outarray[Peta],1)<<" Min ="<<Vmath::Vmin(nSolutionPts, outarray[Peta],1)<<"\t";
         // cout<<"ADP: "<< Vmath::Vmax(nSolutionPts, outarray[ADP],1)<<" Min ="<<Vmath::Vmin(nSolutionPts, outarray[ADP],1)<<"\t";
      
    }
    

/**
 * @brief Write non-field variables PB and SE to files
 * 
 * @param time   Current time.
 */
    void UnsteadyAdvectionDiffusion::OutputPBandSE(NekDouble time)
    {
        ofstream myfile;
        // Number of solution points
        int nSolutionPts = GetNpoints();

        std::map<std::string, int>::const_iterator it;
        for (it = varSENames.begin(); it != varSENames.end(); ++it)
        {   //cout<<"\nWriting SE variable in file";
            std::string filename =  it->first+".csv."+ boost::lexical_cast<std::string>(check);
            //cout<<filename;
            myfile.open(filename.c_str());
            myfile << "x coord;y coord;z coord;scalar\n";
            for(int i = 0; i < nSolutionPts; i++)
                myfile<<coord_0[i]<<";"<<coord_1[i]<<";"<<coord_2[i]<<";"<<se[it->second][i]*NonDse[it->second]<<"\n";//*NonDse[it->second]
            myfile.close();
        }
        for (it = varPBNames.begin(); it != varPBNames.end(); ++it)
        {   
            std::string filename =  it->first+".csv."+ boost::lexical_cast<std::string>(check);
            myfile.open(filename.c_str());
            myfile << "x coord;y coord;z coord;scalar\n";
            for(int i = 0; i < nSolutionPts; i++)
                myfile<<coord_0[i]<<";"<<coord_1[i]<<";"<<coord_2[i]<<";"<<pb[it->second][i]*NonDpb[it->second]<<"\n";//*NonDpb[it->second]
            myfile.close();
        }
        std::string filename = "PLPba.csv."+ boost::lexical_cast<std::string>(check);
        myfile.open(filename.c_str());
        myfile << "x coord;y coord;z coord;scalar\n";
        for(int i = 0; i < nSolutionPts; i++)
            myfile<<coord_0[i]<<";"<<coord_1[i]<<";"<<coord_2[i]<<";"<<PLPba[i]*NonDpb[it->second]<<"\n";//*NonDpb[it->second]
        myfile.close();
        
        filename = "PLPsea.csv."+ boost::lexical_cast<std::string>(check);
        myfile.open(filename.c_str());
        myfile << "x coord;y coord;z coord;scalar\n";
        for(int i = 0; i < nSolutionPts; i++)
            myfile<<coord_0[i]<<";"<<coord_1[i]<<";"<<coord_2[i]<<";"<<PLPsea[i]*NonDpb[it->second]<<"\n";//*NonDpb[it->second]
        myfile.close();

        
        
    }

/**
 * @brief Return the flux vector for the advection part.
 * 
 * @param physfield   Fields.
 * @param flux        Resulting flux.
 */
    void UnsteadyAdvectionDiffusion::GetFluxVectorAdv(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
    {
        ASSERTL1(flux[0].num_elements() == m_velocity.num_elements(),
                 "Dimension of flux array and velocity array do not match");

        const int nSolutionPts = m_fields[0]->GetNpoints();

        for (int i = 0; i < flux.num_elements(); ++i)
        {
            for (int j = 0; j < flux[0].num_elements(); ++j)
            {
                Vmath::Vmul(nSolutionPts, physfield[i], 1, m_velocity[j], 1,
                            flux[i][j], 1);
            }
        }
    }

/**
 * @brief Return the flux vector for the diffusion part.
 *      
 * @param i           Equation number.
 * @param j           Spatial direction.
 * @param physfield   Fields.
 * @param derivatives First order derivatives.
 * @param flux        Resulting flux.
 */
    void UnsteadyAdvectionDiffusion::GetFluxVectorDiff(
        const int i,
        const int j,
        const Array<OneD, Array<OneD, NekDouble> > &physfield,
        Array<OneD, Array<OneD, NekDouble> > &derivatives,
        Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        for (int k = 0; k < flux.num_elements(); ++k)
        {
            Vmath::Zero(GetNpoints(), flux[k], 1);
        }
        Vmath::Vcopy(GetNpoints(), physfield[i], 1, flux[j], 1);
    }
    
    void UnsteadyAdvectionDiffusion::v_GenerateSummary(
        SolverUtils::SummaryList& s)
    {
        AdvectionSystem::v_GenerateSummary(s);
    }

    
/**
 * Perform the extrapolation.
 */
    bool UnsteadyAdvectionDiffusion::v_PreIntegrate(int step)
    {
        if(m_subSteppingScheme)
        {
            SubStepAdvance(m_intSoln,step,m_time);
        }

        return false;
    }


/** 
 * 
 */
    void UnsteadyAdvectionDiffusion::SubStepAdvance(
        const LibUtilities::TimeIntegrationSolutionSharedPtr &integrationSoln, 
        int nstep, 
        NekDouble time)
    {
        int n;
        int nsubsteps;
        
        NekDouble dt; 
        
        Array<OneD, Array<OneD, NekDouble> > fields, velfields;
        
        static int ncalls = 1;
        int  nint         = min(ncalls++, m_intSteps);
        
        Array<OneD, NekDouble> CFL(m_fields[0]->GetExpSize(), 
                                   m_cflSafetyFactor);
        
        LibUtilities::CommSharedPtr comm = m_session->GetComm();

        // Get the proper time step with CFL control
        dt = GetSubstepTimeStep();

        nsubsteps = (m_timestep > dt)? ((int)(m_timestep/dt)+1):1; 
        nsubsteps = max(m_minsubsteps, nsubsteps);

        dt = m_timestep/nsubsteps;
        
        if (m_infosteps && !((nstep+1)%m_infosteps) && comm->GetRank() == 0)
        {
            cout << "Sub-integrating using "<< nsubsteps 
                 << " steps over Dt = "     << m_timestep 
                 << " (SubStep CFL="        << m_cflSafetyFactor << ")"<< endl;
        }

        for (int m = 0; m < nint; ++m)
        {
            // We need to update the fields held by the m_integrationSoln
            fields = integrationSoln->UpdateSolutionVector()[m];
            
            // Initialise NS solver which is set up to use a GLM method
            // with calls to EvaluateAdvection_SetPressureBCs and
            // SolveUnsteadyStokesSystem
            LibUtilities::TimeIntegrationSolutionSharedPtr 
                SubIntegrationSoln = m_subStepIntegrationScheme->
                InitializeScheme(dt, fields, time, m_subStepIntegrationOps);
            
            for(n = 0; n < nsubsteps; ++n)
            {
                fields = m_subStepIntegrationScheme->TimeIntegrate(n, dt, SubIntegrationSoln,
                                                                   m_subStepIntegrationOps);
            }
            
            // Reset time integrated solution in m_integrationSoln 
            integrationSoln->SetSolVector(m,fields);
        }
    }
    

/** 
 * 
 */
    NekDouble UnsteadyAdvectionDiffusion::GetSubstepTimeStep()
    { 
        int n_element = m_fields[0]->GetExpSize(); 

        const Array<OneD, int> ExpOrder=m_fields[0]->EvalBasisNumModesMaxPerExp();
        Array<OneD, int> ExpOrderList (n_element, ExpOrder);
        
        const NekDouble cLambda = 0.2; // Spencer book pag. 317
        
        Array<OneD, NekDouble> tstep      (n_element, 0.0);
        Array<OneD, NekDouble> stdVelocity(n_element, 0.0);

        stdVelocity = GetMaxStdVelocity(m_velocity);
        
        for(int el = 0; el < n_element; ++el)
        {
            tstep[el] = m_cflSafetyFactor / 
                (stdVelocity[el] * cLambda * 
                 (ExpOrder[el]-1) * (ExpOrder[el]-1));
        }
        
        NekDouble TimeStep = Vmath::Vmin(n_element, tstep, 1);
        m_session->GetComm()->AllReduce(TimeStep,LibUtilities::ReduceMin);        
        
        return TimeStep;
    }

    void UnsteadyAdvectionDiffusion::SetUpSubSteppingTimeIntegration(
        int intMethod,
        const LibUtilities::TimeIntegrationWrapperSharedPtr &IntegrationScheme)
    {
        // Set to 1 for first step and it will then be increased in
        // time advance routines
        switch(intMethod)
        {
        case LibUtilities::eBackwardEuler:
        case LibUtilities::eBDFImplicitOrder1: 
        {
            m_subStepIntegrationScheme = LibUtilities::GetTimeIntegrationWrapperFactory().CreateInstance("ForwardEuler");
                
        }
        break;
        case LibUtilities::eBDFImplicitOrder2:
        {
            m_subStepIntegrationScheme = LibUtilities::GetTimeIntegrationWrapperFactory().CreateInstance("RungeKutta2_ImprovedEuler");
        }
        break;
        default:
            ASSERTL0(0,"Integration method not suitable: Options include BackwardEuler or BDFImplicitOrder1");
            break;
        }
        m_intSteps = IntegrationScheme->GetIntegrationSteps();
	
        // set explicit time-integration class operators
        m_subStepIntegrationOps.DefineOdeRhs(&UnsteadyAdvectionDiffusion::SubStepAdvection, this);
        m_subStepIntegrationOps.DefineProjection(&UnsteadyAdvectionDiffusion::SubStepProjection, this);
    }
    
/** 
 * Explicit Advection terms used by SubStepAdvance time integration
 */
    void UnsteadyAdvectionDiffusion::SubStepAdvection(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,  
        Array<OneD, Array<OneD,       NekDouble> > &outarray,
        const NekDouble time)
    {
        int i;
        int nVariables     = inarray.num_elements();
        
        /// Get the number of coefficients
        int ncoeffs = m_fields[0]->GetNcoeffs(); 
        
        /// Define an auxiliary variable to compute the RHS 
        Array<OneD, Array<OneD, NekDouble> > WeakAdv(nVariables);
        WeakAdv[0] = Array<OneD, NekDouble> (ncoeffs*nVariables);
        for(i = 1; i < nVariables; ++i)
        {
            WeakAdv[i] = WeakAdv[i-1] + ncoeffs;
        }
        
        // Currently assume velocity field is time independent and does not therefore
        // need extrapolating. 
        // RHS computation using the advection base class
        m_advObject->Advect(nVariables, m_fields, m_velocity, 
                            inarray, outarray, time);

        for(i = 0; i < nVariables; ++i)
        {
            m_fields[i]->IProductWRTBase(outarray[i],WeakAdv[i]); 
            // negation requried due to sign of DoAdvection term to be consistent
            Vmath::Neg(ncoeffs, WeakAdv[i], 1);
        }
        
        AddAdvectionPenaltyFlux(m_velocity, inarray, WeakAdv);

        
        /// Operations to compute the RHS
        for(i = 0; i < nVariables; ++i)
        {
            // Negate the RHS
            Vmath::Neg(ncoeffs, WeakAdv[i], 1);

            /// Multiply the flux by the inverse of the mass matrix
            m_fields[i]->MultiplyByElmtInvMass(WeakAdv[i], WeakAdv[i]);
            
            /// Store in outarray the physical values of the RHS
            m_fields[i]->BwdTrans(WeakAdv[i], outarray[i]);
        }
    }
        
/** 
 * Projection used by SubStepAdvance time integration
 */
    void UnsteadyAdvectionDiffusion::SubStepProjection(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,  
        Array<OneD, Array<OneD, NekDouble> > &outarray, 
        const NekDouble time)
    {
        ASSERTL1(inarray.num_elements() == outarray.num_elements(),"Inarray and outarray of different sizes ");

        for(int i = 0; i < inarray.num_elements(); ++i)
        {
            Vmath::Vcopy(inarray[i].num_elements(),inarray[i],1,outarray[i],1);
        }
    }

    void UnsteadyAdvectionDiffusion::AddAdvectionPenaltyFlux(
        const Array<OneD, const Array<OneD, NekDouble> > &velfield, 
        const Array<OneD, const Array<OneD, NekDouble> > &physfield, 
        Array<OneD, Array<OneD, NekDouble> > &Outarray)
    {
        ASSERTL1(physfield.num_elements() == Outarray.num_elements(),
                 "Physfield and outarray are of different dimensions");
        
        int i;
        
        /// Number of trace points
        int nTracePts   = m_fields[0]->GetTrace()->GetNpoints();

        /// Forward state array
        Array<OneD, NekDouble> Fwd(3*nTracePts);
        
        /// Backward state array
        Array<OneD, NekDouble> Bwd = Fwd + nTracePts;

        /// upwind numerical flux state array
        Array<OneD, NekDouble> numflux = Bwd + nTracePts;
        
        /// Normal velocity array
        Array<OneD, NekDouble> Vn  = GetNormalVel(velfield);
        
        for(i = 0; i < physfield.num_elements(); ++i)
        {
            /// Extract forwards/backwards trace spaces
            /// Note: Needs to have correct i value to get boundary conditions
            m_fields[i]->GetFwdBwdTracePhys(physfield[i], Fwd, Bwd);
            
            /// Upwind between elements
            m_fields[0]->GetTrace()->Upwind(Vn, Fwd, Bwd, numflux);

            /// Construct difference between numflux and Fwd,Bwd
            Vmath::Vsub(nTracePts, numflux, 1, Fwd, 1, Fwd, 1);
            Vmath::Vsub(nTracePts, numflux, 1, Bwd, 1, Bwd, 1);

            /// Calculate the numerical fluxes multipling Fwd, Bwd and
            /// numflux by the normal advection velocity
            Vmath::Vmul(nTracePts, Fwd, 1, Vn, 1, Fwd, 1);
            Vmath::Vmul(nTracePts, Bwd, 1, Vn, 1, Bwd, 1);

            m_fields[0]->AddFwdBwdTraceIntegral(Fwd,Bwd,Outarray[i]);
        }
    }


    Array<OneD, NekDouble> UnsteadyAdvectionDiffusion::GetMaxStdVelocity(
        const Array<OneD, Array<OneD,NekDouble> > inarray)
    {
        
        int n_points_0      = m_fields[0]->GetExp(0)->GetTotPoints();
        int n_element       = m_fields[0]->GetExpSize();       
        int nvel            = inarray.num_elements();
        int cnt; 

        ASSERTL0(nvel >= 2, "Method not implemented for 1D");
        
        NekDouble pntVelocity;
        
        // Getting the standard velocity vector on the 2D normal space
        Array<OneD, Array<OneD, NekDouble> > stdVelocity(nvel);
        Array<OneD, NekDouble> maxV(n_element, 0.0);
        LibUtilities::PointsKeyVector ptsKeys;
        
        for (int i = 0; i < nvel; ++i)
        {
            stdVelocity[i] = Array<OneD, NekDouble>(n_points_0);
        }
        
        if (nvel == 2)
        {
            cnt = 0.0;
            for (int el = 0; el < n_element; ++el)
            { 
                int n_points = m_fields[0]->GetExp(el)->GetTotPoints();
                ptsKeys = m_fields[0]->GetExp(el)->GetPointsKeys();
                
                // reset local space if necessary
                if(n_points != n_points_0)
                {
                    for (int j = 0; j < nvel; ++j)
                    {
                        stdVelocity[j] = Array<OneD, NekDouble>(n_points);
                    }
                    n_points_0 = n_points;
                }		
                
                Array<TwoD, const NekDouble> gmat = 
                    m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()->GetDerivFactors(ptsKeys);
                
                if (m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()->GetGtype()
                    == SpatialDomains::eDeformed)
                {
                    for (int i = 0; i < n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][i]*inarray[0][i+cnt] 
                            + gmat[2][i]*inarray[1][i+cnt];
                        
                        stdVelocity[1][i] = gmat[1][i]*inarray[0][i+cnt] 
                            + gmat[3][i]*inarray[1][i+cnt];
                    }
                }
                else
                {
                    for (int i = 0; i < n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][0]*inarray[0][i+cnt] 
                            + gmat[2][0]*inarray[1][i+cnt];
                        
                        stdVelocity[1][i] = gmat[1][0]*inarray[0][i+cnt] 
                            + gmat[3][0]*inarray[1][i+cnt];
                    }
                }
                
                cnt += n_points;
                
                
                for (int i = 0; i < n_points; i++)
                {
                    pntVelocity = stdVelocity[0][i]*stdVelocity[0][i] 
                        + stdVelocity[1][i]*stdVelocity[1][i];
                    
                    if (pntVelocity>maxV[el])
                    {
                        maxV[el] = pntVelocity;
                    }
                }
                maxV[el] = sqrt(maxV[el]);
            }
        }
        else
        {
            cnt = 0;
            for (int el = 0; el < n_element; ++el)
            { 
                
                int n_points = m_fields[0]->GetExp(el)->GetTotPoints();
                ptsKeys = m_fields[0]->GetExp(el)->GetPointsKeys();
                
                // reset local space if necessary
                if(n_points != n_points_0)
                {
                    for (int j = 0; j < nvel; ++j)
                    {
                        stdVelocity[j] = Array<OneD, NekDouble>(n_points);
                    }
                    n_points_0 = n_points;
                }		
                
                Array<TwoD, const NekDouble> gmat =
                    m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()->GetDerivFactors(ptsKeys);
                
                if (m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()->GetGtype()
                    == SpatialDomains::eDeformed)
                {
                    for (int i = 0; i < n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][i]*inarray[0][i+cnt] 
                            + gmat[3][i]*inarray[1][i+cnt] 
                            + gmat[6][i]*inarray[2][i+cnt];
                        
                        stdVelocity[1][i] = gmat[1][i]*inarray[0][i+cnt] 
                            + gmat[4][i]*inarray[1][i+cnt] 
                            + gmat[7][i]*inarray[2][i+cnt];
                        
                        stdVelocity[2][i] = gmat[2][i]*inarray[0][i+cnt] 
                            + gmat[5][i]*inarray[1][i+cnt] 
                            + gmat[8][i]*inarray[2][i+cnt];
                    }
                }
                else
                {
                    for (int i = 0; i < n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][0]*inarray[0][i+cnt] 
                            + gmat[3][0]*inarray[1][i+cnt] 
                            + gmat[6][0]*inarray[2][i+cnt];
                        
                        stdVelocity[1][i] = gmat[1][0]*inarray[0][i+cnt] 
                            + gmat[4][0]*inarray[1][i+cnt] 
                            + gmat[7][0]*inarray[2][i+cnt];
                        
                        stdVelocity[2][i] = gmat[2][0]*inarray[0][i+cnt] 
                            + gmat[5][0]*inarray[1][i+cnt] 
                            + gmat[8][0]*inarray[2][i+cnt];
                    }
                }
                
                cnt += n_points;
                
                for (int i = 0; i < n_points; i++)
                {
                    pntVelocity = stdVelocity[0][i]*stdVelocity[0][i] 
                        + stdVelocity[1][i]*stdVelocity[1][i] 
                        + stdVelocity[2][i]*stdVelocity[2][i];
                    
                    if (pntVelocity > maxV[el])
                    {
                        maxV[el] = pntVelocity;
                    }
                }

                maxV[el] = sqrt(maxV[el]);
                //couty << maxV[el]*maxV[el] << endl;
            }
        }
		
        return maxV;
    }
}
