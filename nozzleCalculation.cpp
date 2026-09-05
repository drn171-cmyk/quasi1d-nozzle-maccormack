// This code for solving differantial equations of quasi 1D flow through nozzle by using MacCormak method


//Importing libraries

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

double timeStepSize(double u, double T, double dx) {
    
    double courantNumber = 0.1; // Lowered from 0.5: the shocked case has a much
    // steeper initial transient (back pressure imposed on an initial guess built
    // for the fully-supersonic case), so a smaller CFL number is needed for
    // stability. 0.5 is fine again if you go back to the shockless case.
    
    double a = sqrt(T); // Non dimentioalised speed of sound
    
    return courantNumber * dx / (u + a);
}

int main() {

    // Constants

    double gamma = 1.4; // The ratio of spefic heat
    double Cp = 1005; // Specific heat of air
    double R = 287; // Gas constant of air

    // Back pressure (non-dimensionalised by reservoir pressure p0) imposed at the
    // outlet. This value corresponds to Anderson's classic case with a normal
    // shock standing inside the diverging section (subsonic exit). Set to a
    // value <= the fully-supersonic exit pressure (~0.0225 at N=1000) to
    // recover the original shockless/fully-supersonic case instead.
    double pBack = 0.6784;

    // Artificial viscosity coefficient. The shocked/subsonic-exit case has a
    // discontinuity (the shock) in the interior of the domain, which the
    // basic MacCormack scheme cannot capture cleanly without added damping:
    // without this the scheme blows up (oscillations grow until values go
    // negative/NaN near the shock). Typical value from Anderson's book.
    double Cx = 0.3;
    
    // Defining grid

    int N = 1000; // Grid number    
    double L = 3; // Nozzle length
    int maxTimeStepNumber = 400000; // Safety cap; loop normally stops earlier via residual check
    double convergenceTolerance = 1e-8; // Stop when max relative change in rho drops below this
    double dx = L/(N-1); // Length of the each grid
    
    // Defining Vectors

    vector<double> x(N); // Position vector
    vector<double> A(N); // Cross section vector
    vector<double> rho(N); // Density vector
    vector<double> u(N); // Velocity vector
    vector<double> p(N); // Static pressure vector
    vector<double> T(N); // Static temperature vector
    vector<double> U1(N), U2(N), U3(N); // U is the solution vector
    vector<double> F1(N), F2(N), F3(N); // F is the flux vector
    vector<double> S2(N); // S is the source vector. S1 and S3 equal 0. Therefore, only S2 is defined
    vector<double> U1Star(N), U2Star(N), U3Star(N); // Prediction vectors for U
    vector<double> F1Star(N), F2Star(N), F3Star(N); // Prediction vectors for F
    vector<double> S2Star(N); // Prediction vectors for F
    vector<double> pStarVec(N); // Predicted pressure, kept for the artificial viscosity stencil

    // Loop for defining initial conditions

    for(int i = 0; i < N; i++){

        x[i] = i * dx; // Position
        A[i] = 1.0 + 2.2 * pow((x[i] - 1.5), 2); // Cross section

        // This initial conditions are taken from Anderson's Modern Compressible Flow
        // p, T, and rho is non dimentioalised.
        
        rho[i] = 1.0 - 0.3146 * x[i]; // Density
        T[i] = 1.0 - 0.2314 * x[i]; // Static Temperature
        u[i] = ( 0.1 + 1.09 * x[i]) * sqrt(T[i]); // Velocity 

        // Initialsing U vectors

        U1[i] = rho[i] * A[i]; // U1 = rho * A
        U2[i] = u[i] * rho[i] * A[i]; // U2 = u * rho * A
        U3[i] = rho[i] * (T[i] / (gamma - 1.0) + (gamma / 2.0) * pow(u[i], 2)) * A[i]; // U3 = rho * (Internal Energy + Kinetic Energy) * A
    }

    int t;
    for( t = 0; t < maxTimeStepNumber; t++ ){

        // dt must be the MINIMUM over the whole grid, not just the outlet.
        // With a fully supersonic exit, velocity increases monotonically so
        // the outlet was always the most restrictive point. Once the exit
        // is forced subsonic (fixed back pressure), that's no longer true:
        // the fastest point is now somewhere inside the nozzle, and using
        // only the outlet gives a dt that's too large there, which is what
        // was causing the blow-up near the shock/exit.
        double dt = timeStepSize(u[0], T[0], dx);
        for (int i = 1; i < N; i++) {
            dt = min(dt, timeStepSize(u[i], T[i], dx));
        }

        vector<double> rhoOld = rho; // Kept to measure convergence (residual) after this step

        // Predictor step

        // Calculate F and S

        for (int i = 0; i < N; i++) {

        F1[i] = U2[i];
        F2[i] = (pow(U2[i], 2) / U1[i]) + ((gamma - 1.0)/gamma) * (U3[i] - (gamma/2.0)*pow(U2[i], 2)/U1[i]);
        F3[i] = gamma * (U2[i] * U3[i] / U1[i]) - (gamma * (gamma - 1.0) / 2.0) * (pow(U2[i], 3) / pow(U1[i], 2));

        p[i] = ((gamma - 1.0)/gamma) * (U3[i] - (gamma/2.0) * pow(U2[i], 2) / U1[i]) / A[i]; // p = rho*T

        // S2 = p * dA/dx, dA/dx is calculated via center difference

        if (i > 0 && i < N - 1) {
           
            S2[i] = p[i] * (A[i+1] - A[i-1]) / (2.0 * dx);
        
        }

        }

        S2[0] = 0.0; S2[N-1] = 0.0; // Source term is equal zero at boundary points.

        // U_Star calculation
        
        for (int i = 1; i < N - 1; i++) {
            
            U1Star[i] = U1[i] - (dt / dx) * (F1[i+1] - F1[i]);
            U2Star[i] = U2[i] - (dt / dx) * (F2[i+1] - F2[i]) + dt * S2[i];
            
            U3Star[i] = U3[i] - (dt / dx) * (F3[i+1] - F3[i]);
        }

        U1Star[0] = U1[0]; U2Star[0] = U2[0]; U3Star[0] = U3[0];
        U1Star[N-1] = U1[N-1]; U2Star[N-1] = U2[N-1]; U3Star[N-1] = U3[N-1];

        // Artificial viscosity (predictor step), using pressure at time
        // level n. Damps the odd/even oscillations that otherwise grow
        // without bound around a captured shock.
        for (int i = 1; i < N - 1; i++) {
            double denom = p[i+1] + 2.0 * p[i] + p[i-1];
            double visc = Cx * fabs(p[i+1] - 2.0 * p[i] + p[i-1]) / denom;
            U1Star[i] += visc * (U1[i+1] - 2.0 * U1[i] + U1[i-1]);
            U2Star[i] += visc * (U2[i+1] - 2.0 * U2[i] + U2[i-1]);
            U3Star[i] += visc * (U3[i+1] - 2.0 * U3[i] + U3[i-1]);
        }

        // Corrector Step

        // Calculate F_star ve S_star by using U_Star
        
        for (int i = 0; i < N; i++) {
            
            F1Star[i] = U2Star[i];
            F2Star[i] = (pow(U2Star[i], 2) / U1Star[i]) + ((gamma - 1.0) / gamma) * (U3Star[i] - (gamma / 2.0) * pow(U2Star[i], 2) / U1Star[i]);
            F3Star[i] = gamma * (U2Star[i] * U3Star[i] / U1Star[i]) - (gamma * (gamma - 1.0) / 2.0) * (pow(U2Star[i], 3) / pow(U1Star[i], 2));

            double pStar = ((gamma - 1.0) / gamma) * (U3Star[i] - (gamma / 2.0) * pow(U2Star[i], 2) / U1Star[i]) / A[i];
            pStarVec[i] = pStar; // Saved for the artificial viscosity stencil below (scale is irrelevant there since it cancels in the |.|/(.) ratio)
            
            // Backward difference for S2Star 
            if (i > 0) {
                S2Star[i] = pStar * (A[i] - A[i-1]) / dx;
            }
        }

        // Update U values

        for (int i = 1; i < N - 1; i++) {
            
            double dU1_dt_corr = -(F1Star[i] - F1Star[i-1]) / dx;
            double dU2_dt_corr = -(F2Star[i] - F2Star[i-1]) / dx + S2Star[i];
            double dU3_dt_corr = -(F3Star[i] - F3Star[i-1]) / dx;

            U1[i] = 0.5 * (U1[i] + U1Star[i] + dU1_dt_corr * dt);
            U2[i] = 0.5 * (U2[i] + U2Star[i] + dU2_dt_corr * dt);
            U3[i] = 0.5 * (U3[i] + U3Star[i] + dU3_dt_corr * dt);
        }

        // Artificial viscosity (corrector step), using the predicted
        // pressure pStarVec. Same purpose as the predictor-step version
        // above, applied to the corrected values.
        for (int i = 1; i < N - 1; i++) {
            double denom = pStarVec[i+1] + 2.0 * pStarVec[i] + pStarVec[i-1];
            double visc = Cx * fabs(pStarVec[i+1] - 2.0 * pStarVec[i] + pStarVec[i-1]) / denom;
            U1[i] += visc * (U1Star[i+1] - 2.0 * U1Star[i] + U1Star[i-1]);
            U2[i] += visc * (U2Star[i+1] - 2.0 * U2Star[i] + U2Star[i-1]);
            U3[i] += visc * (U3Star[i+1] - 2.0 * U3Star[i] + U3Star[i-1]);
        }

        // Boundary Conditions

        // Inlet 
        // Due to non dimentionalsied reservouir temperature and denstiy, they are equal to 1.
        
        U1[0] = A[0]; 
        U2[0] = 2.0 * U2[1] - U2[2]; // Ekstrapolasion from forward nodes
        double u_inlet = U2[0] / U1[0];
        U3[0] = U1[0] * (1.0 / (gamma - 1.0) + (gamma / 2.0) * pow(u_inlet, 2));

        // Outlet
        // Shocked / subsonic-exit case: a normal shock forms inside the
        // diverging section, so the exit flow is subsonic and an acoustic
        // signal can travel upstream from the ambient into the nozzle.
        // Only rho and u are extrapolated from the interior; pressure at
        // the exit is fixed to the imposed back pressure pBack, and T, then
        // U1/U2/U3, are reconstructed from these primitive variables.

        double rho_im1 = U1[N-2] / A[N-2];
        double rho_im2 = U1[N-3] / A[N-3];
        double u_im1   = U2[N-2] / U1[N-2];
        double u_im2   = U2[N-3] / U1[N-3];

        double rho_exit = 2.0 * rho_im1 - rho_im2; // Extrapolated
        double u_exit   = 2.0 * u_im1   - u_im2;   // Extrapolated
        double p_exit   = pBack;                   // Fixed back pressure
        double T_exit   = p_exit / rho_exit;        // From p = rho * T

        U1[N-1] = rho_exit * A[N-1];
        U2[N-1] = rho_exit * u_exit * A[N-1];
        U3[N-1] = rho_exit * (T_exit / (gamma - 1.0) + (gamma / 2.0) * pow(u_exit, 2)) * A[N-1];

        
        // Update rho, u, T, and p
        
        for (int i = 0; i < N; i++) {
            rho[i] = U1[i] / A[i];
            u[i]   = U2[i] / U1[i];
            T[i]   = (gamma - 1.0) * (U3[i] / U1[i] - (gamma / 2.0) * pow(u[i], 2));
            p[i]   = rho[i] * T[i];
        }

        // Convergence check: steady state is reached once d(rho)/dt is
        // negligible everywhere. This makes the stopping point independent
        // of the grid size N (dt shrinks as N grows, so a fixed iteration
        // count that works for a coarse grid stops far too early on a fine
        // one, as happened here with N = 1000).

        double maxDrhoDt = 0.0;
        for (int i = 0; i < N; i++) {
            double drhodt = fabs(rho[i] - rhoOld[i]) / dt;
            if (drhodt > maxDrhoDt) maxDrhoDt = drhodt;
        }

        if (t % 100 == 0) {
            cout << "Iteration " << t << ", dt = " << dt
                 << ", max d(rho)/dt = " << maxDrhoDt << endl;
        }

        if (maxDrhoDt < convergenceTolerance) {
            cout << "Converged after " << t << " iterations "
                 << "(max d(rho)/dt = " << maxDrhoDt << ")." << endl;
            break;
        }
    }

    if (t >= maxTimeStepNumber) {
        cout << "Warning: reached the " << maxTimeStepNumber
             << " iteration cap without fully converging. "
             << "Consider raising maxTimeStepNumber or loosening convergenceTolerance."
             << endl;
    }


    // Writing results on a CVS file
    
    ofstream outFile("nozzleResults.csv"); // Create file
    
    if (outFile.is_open()) {
        
        // Writes column head

        outFile << "x,A,rho,u,T,p,Mach\n";
        
        // Writes final values of all grid

        for (int i = 0; i < N; i++) {
            double Mach = u[i] / sqrt(T[i]); // Calculate and add Mach number
            outFile << x[i] << "," 
                    << A[i] << "," 
                    << rho[i] << "," 
                    << u[i] << "," 
                    << T[i] << "," 
                    << p[i] << "," 
                    << Mach << "\n";
        }
        
        outFile.close();
        cout << "Calculation complated. Results are saved to 'nozzleResults.csv' file." << endl;
    
    } else {
        
        cout << "File cannot created" << endl;
    
    }


    return 0;
}