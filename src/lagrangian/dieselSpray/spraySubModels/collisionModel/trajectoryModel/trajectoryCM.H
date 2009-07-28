vector v1 = p1().U();
vector v2 = p2().U();
scalar pc = spray_.p()[p1().cell()];

vector vRel = v1 - v2;
scalar magVRel = mag(vRel);

vector p = p2().position() - p1().position();
scalar dist = mag(p);

scalar vAlign = vRel & (p/(dist+SMALL));

if (vAlign > 0)
{
    scalar sumD = p1().d() + p2().d();

    if (vAlign*dt > dist - 0.5*sumD)
    {
        scalar v1Mag = mag(v1);
        scalar v2Mag = mag(v2);
        vector nv1 = v1/v1Mag;
        vector nv2 = v2/v2Mag;

        scalar v1v2 = nv1 & nv2;
        scalar v1p = nv1 & p;
        scalar v2p = nv2 & p;

        scalar det = 1.0 - v1v2*v1v2;

        scalar alpha = 1.0e+20;
        scalar beta = 1.0e+20;

        if (mag(det) > 1.0e-4)
        {
            beta = -(v2p - v1v2*v1p)/det;
            alpha = v1p + v1v2*beta;
        }

        alpha /= v1Mag*dt;
        beta /= v2Mag*dt;

        // is collision possible within this timestep
        if ((alpha>0) && (alpha<1.0) && (beta>0) && (beta<1.0)) 
        {
            vector p1c = p1().position() + alpha*v1*dt;
            vector p2c = p2().position() + beta*v2*dt;

            scalar closestDist = mag(p1c-p2c);

            scalar collProb = 
                pow(0.5*sumD/max(0.5*sumD, closestDist), cSpace_)
               *exp(-cTime_*mag(alpha-beta));

            scalar xx = rndGen_.scalar01();

            spray::iterator pMin = p1;
            spray::iterator pMax = p2;

            scalar dMin = pMin().d();
            scalar dMax = pMax().d();

            if (dMin > dMax)
            {
                dMin = pMax().d();
                dMax = pMin().d();
                pMin = p2;
                pMax = p1;
            }

            scalar rhoMax = spray_.fuels().rho(pc, pMax().T(), pMax().X());
            scalar rhoMin = spray_.fuels().rho(pc, pMin().T(), pMin().X());
            scalar mMax = pMax().m();
            scalar mMin = pMin().m();
            scalar nMax = pMax().N(rhoMax);
            scalar nMin = pMin().N(rhoMin);

            scalar mdMin = mMin/nMin;

            // collision occur
            if ((xx < collProb) && (mMin > VSMALL) && (mMax > VSMALL))
            {
                scalar mTot = mMax + mMin;

                scalar gamma = dMax/max(dMin, 1.0e-12);
                scalar f = gamma*gamma*gamma + 2.7*gamma - 2.4*gamma*gamma;

                vector momMax = mMax*pMax().U();
                vector momMin = mMin*pMin().U();

                // use mass-averaged temperature to calculate We number
                scalar averageTemp = (pMax().T()*mMax + pMin().T()*mMin)/mTot;
                // and mass averaged mole fractions ...
                scalarField 
                    Xav((pMax().m()*pMax().X()+pMin().m()*pMin().X())
                   /(pMax().m() + pMin().m()));

                scalar sigma = spray_.fuels().sigma(pc, averageTemp, Xav);
                sigma = max(1.0e-6, sigma);
                scalar rho = spray_.fuels().rho(pc, averageTemp, Xav);

                scalar dMean = sqrt(dMin*dMax);
                scalar WeColl = 
                    max(1.0e-12, 0.5*rho*magVRel*magVRel*dMean/sigma);

                // coalescence only possible when parcels are close enoug

                scalar coalesceProb = min(1.0, 2.4*f/WeColl);

                scalar prob = rndGen_.scalar01();

                // Coalescence
                if ( prob < coalesceProb && coalescence_) 
                {
                    // How 'many' of the droplets coalesce
                    scalar nProb = prob*nMin/nMax;

                    // Conservation of mass, momentum and energy

                    pMin().m() -= nMax*nProb*mdMin;

                    scalar newMinMass = pMin().m();
                    scalar newMaxMass = mMax + (mMin - newMinMass);
                    pMax().m() = newMaxMass;

                    pMax().T() = 
                        (averageTemp*mTot - newMinMass*pMin().T())/newMaxMass;
                    rhoMax = spray_.fuels().rho(pc, pMax().T(), pMax().X());

                    pMax().d() = 
                        pow
                        (
                            6.0*newMaxMass/(rhoMax*mathematicalConstant::pi*nMax),
                            1.0/3.0
                        );

                    pMax().U() = 
                        (momMax + (1.0-newMinMass/mMin)*momMin)/newMaxMass;

                    // update the liquid molar fractions
                    scalarField Ymin = spray_.fuels().Y(pMin().X());
                    scalarField Ymax = spray_.fuels().Y(pMax().X());
                    scalarField Ynew = mMax*Ymax + (mMin - newMinMass)*Ymin;
                    scalar Wlinv = 0.0;
                    forAll(Ynew, i)
                    {
                        Wlinv += Ynew[i]/spray_.fuels().properties()[i].W();
                    }
                    forAll(Ynew, i)
                    {
                        pMax().X()[i] = 
                            Ynew[i]/(spray_.fuels().properties()[i].W()*Wlinv);
                    }


                }
                // Grazing collision (no coalescence)
                else
                {
                    scalar gf = sqrt(prob) - sqrt(coalesceProb);
                    scalar denom = 1.0 - sqrt(coalesceProb);
                    if (denom < 1.0e-5) {
                        denom = 1.0;
                    }
                    gf /= denom;

                    // if gf negative, this means that coalescence is turned off
                    // and these parcels should have coalesced
                    gf = max(0.0, gf);

                    scalar rho1 = spray_.fuels().rho(pc, p1().T(), p1().X());
                    scalar rho2 = spray_.fuels().rho(0.0, p2().T(), p2().X());
                    scalar m1 = p1().m();
                    scalar m2 = p2().m();
                    scalar n1 = p1().N(rho1);
                    scalar n2 = p2().N(rho2);

                    // gf -> 1 => v1p -> p1().U() ...
                    // gf -> 0 => v1p -> momentum/(m1+m2)

                    vector mr = m1*v1 + m2*v2;
                    vector v1p = (mr + m2*gf*vRel)/(m1+m2);
                    vector v2p = (mr - m1*gf*vRel)/(m1+m2);

                    if (n1 < n2)
                    {
                        p1().U() = v1p;
                        p2().U() = (n1*v2p + (n2-n1)*v2)/n2;
                    }
                    else
                    {
                        p1().U() = (n2*v1p + (n1-n2)*v1)/n1;
                        p2().U() = v2p;
                    }

                } // if - coalescence or not

            } // if - collision

        } // if - possible collision (alpha, beta) in timeinterval

    } // if - travelled distance is larger distance between parcels

} // if - parcel travel towards eachother
