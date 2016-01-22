#ifndef PARTICLEFNS_H
#define PARTICLEFNS_H

#include <math.h>

// define global vars
double capillary_center[3] = {0, GRID_LENGTH/2., GRID_LENGTH/2.};

// qsort compare function for sorting by axial distance from capillary
int axial_dist_compare(const void *p1, const void *p2)
{
    const Particle *p_1 = (Particle*)p1;
    const Particle *p_2 = (Particle*)p2;

    double tx = p_1->x - capillary_center[0];
    double ty = p_1->y - capillary_center[1];
    double tz = p_1->z - capillary_center[2];

    // first position vector squared
    double r1 = tx*tx + ty*ty + tz*tz;

    tx = p_2->x - capillary_center[0];
    ty = p_2->y - capillary_center[1];
    tz = p_2->z - capillary_center[2];
    double r2 = tx*tx + ty*ty + tz*tz;

    double res = r1-r2;
    if(fabs(res) < 1e-6)
        return 0;
    else if(res > 1e-6)
        return 1;
    else
        return -1;
}

bool isParticleInDomain(const Particle p)
{
    // assuming that domain lengths can never go -ve
    // or beyond the max GRID_LENGTH
    return ( (p.x >= 0) && (p.x <= GRID_LENGTH) )
        && ( (p.y >= 0) && (p.y <= GRID_LENGTH) )
        && ( (p.z >= 0) && (p.z <= GRID_LENGTH) );
}

Particle randomizeParticleAttribs(Particle inputParticle)
{
    double randNum = (rand() / (double)(RAND_MAX)) * 2. * M_PI;

    double temp1 = (inputParticle.y - (GRID_LENGTH/2.));
    double temp2 = (inputParticle.z - (GRID_LENGTH/2.));

    double rYZ = sqrt(temp1*temp1 + temp2*temp2);

    // adjust for the fact that our origin is at corner
    // whereas MD data origin is at capillary center
    inputParticle.y = rYZ * cos(randNum) + (GRID_LENGTH/2.);
    inputParticle.z = rYZ * sin(randNum) + (GRID_LENGTH/2.);

    // choose a different random number for velocity - just for more randomization
    randNum = (rand() / (double)(RAND_MAX)) * 2. * M_PI;
    double velYZ = sqrt(inputParticle.Vy*inputParticle.Vy + inputParticle.Vz*inputParticle.Vz);
    inputParticle.Vy = velYZ * cos(randNum);
    inputParticle.Vz = velYZ * sin(randNum);

    return inputParticle;
}

void getCellWeights(const double x, const double y, const double z,
                    const double spacing, const double invSpacing,
                    int indices[][3], double *weightFactors)
{
        // get the nearest indices - these must be the lower
        // corner indices of a cell - due to nature of (int) cast
        int iPos = (int)(x*invSpacing);
        int jPos = (int)(y*invSpacing);
        int kPos = (int)(z*invSpacing);

        double hx = (x - iPos*spacing)*invSpacing;
        double hy = (y - jPos*spacing)*invSpacing;
        double hz = (z - kPos*spacing)*invSpacing;

        // ORDERING of filling in the weightfactors and indices
        // X, Y, Z   --> X+1, Y, Z   --> X, Y+1, Z   --> X+1, Y+1, Z
        // X, Y, Z+1 --> X+1, Y, Z+1 --> X, Y+1, Z+1 --> X+1, Y+1, Z+1

        // manually fill in the distances and indices to each of the vertices
        int ti = iPos, tj = jPos, tk = kPos;
        indices[0][0] = ti; indices[0][1] = tj; indices[0][2] = tk;
        weightFactors[0] = (1-hx)*(1-hy)*(1-hz);

        ti = iPos+1, tj = jPos, tk = kPos;
        indices[1][0] = ti; indices[1][1] = tj; indices[1][2] = tk;
        weightFactors[1] = (hx)*(1-hy)*(1-hz);

        ti = iPos; tj = jPos+1; tk = kPos;
        indices[2][0] = ti; indices[2][1] = tj; indices[2][2] = tk;
        weightFactors[2] = (1-hx)*(hy)*(1-hz);

        ti = iPos+1; tj = jPos+1; tk = kPos;
        indices[3][0] = ti; indices[3][1] = tj; indices[3][2] = tk;
        weightFactors[3] = (hx)*(hy)*(1-hz);

        /*********************************************/
        /*********************************************/
        ti = iPos; tj = jPos; tk = kPos+1;
        indices[4][0] = ti; indices[4][1] = tj; indices[4][2] = tk;
        weightFactors[4] = (1-hx)*(1-hy)*(hz);

        ti = iPos+1; tj = jPos; tk = kPos+1;
        indices[5][0] = ti; indices[5][1] = tj; indices[5][2] = tk;
        weightFactors[5] = (hx)*(1-hy)*(hz);

        ti = iPos; tj = jPos+1; tk = kPos+1;
        indices[6][0] = ti; indices[6][1] = tj; indices[6][2] = tk;
        weightFactors[6] = (1-hx)*(hy)*(hz);

        ti = iPos+1; tj = jPos+1; tk = kPos+1;
        indices[7][0] = ti; indices[7][1] = tj; indices[7][2] = tk;
        weightFactors[7] = (hx)*(hy)*(hz);
}

void updateChargeFractions(
         const Particle *particleList, int particleCount,
         double *chargeFractions,
         GridInfo *gInfo
         )
{
    int i, validNodesNumber = 0;
    int searchIndex = 0;

    int numNodes = gInfo->numNodes;
    double invSpacing = gInfo->invSpacing;
    double spacing = gInfo->spacing;

    // arrays for storing weightFactors and indices
    double weightFactors[8];
    int indices[8][3];

    // cache the multiplication factor
    // taking spacing factors to rhs and taking vol as spacing^3
    const double invCellVol = invSpacing*invSpacing*invSpacing;
    const double multFactor = -invCellVol * ELECTRONIC_CHARGE / FREE_SPACE_PERMITTIVITY;

    for(i = 0; i < particleCount; i++)
    {
        const Particle *p = &particleList[i];

        // get the relevant weight factors
        getCellWeights(p->x, p->y, p->z,
                       spacing, invSpacing,
                       indices, weightFactors);

        // update the charge fractions, based on
        // the weight factors
        int c;
        for(c = 0; c < 8; c++)
        {
            int ti = indices[c][0], tj = indices[c][1], tk = indices[c][2];

            // NECESSARY?
            // if within interior points of the grid
            // since we don't want to overwrite boundary
            // portion of the grid
            if( ((ti > 0) && (ti < numNodes-1))
                          &&
                ((tj > 0) && (tj < numNodes-1))
                          &&
                ((tk > 0) && (tk < numNodes-1)) )
            {
                int index = INDEX_1D(numNodes, ti, tj, tk);
                chargeFractions[ index ] += p->charge * weightFactors[c] * multFactor;
            }
        }

        // now loop through the weight factors
        // and update positions in the charge fractions list
        /*
        int j;
        for(j = 0; j < 8; j++)
        {
            ti = indices[j][0];
            tj = indices[j][1];
            tk = indices[j][2];

            // if within interior points of the grid
            // since we don't want to overwrite boundary
            // portion of the grid
            if( ((ti > 0) && (ti < numNodes-1))
                          &&
                ((tj > 0) && (tj < numNodes-1))
                          &&
                ((tk > 0) && (tk < numNodes-1)) )
            {
                int pos = INDEX_1D(numNodes, ti, tj, tk);

                int k;
                for(k = 0; k < validNodesNumber; k++)
                {
                    // wrap around index
                    // reusing last searched value to speed things up?
                    searchIndex = searchIndex % validNodesNumber;
                    if(indexVals[searchIndex].index == pos)
                    {
                        indexVals[searchIndex].val += p->charge * multFactor;
                        break;
                    }
                    searchIndex++;
                }

                // if pos was not found, then write into and expand size
                if(k == validNodesNumber)
                {
                    indexVals[validNodesNumber].index = pos;
                    indexVals[validNodesNumber].val = p->charge * multFactor;
                    validNodesNumber++;
                }
            } // end of if check for indices within valid limits
        } // end of for loop for min distance
        */

    } // end of loop through particles array
}

/*!
 * This routine releases the a specified number of particles into the domain. It uses the lostParticles index
 * count array to keep track of indices of particles which have left the domain. The outline is as follows:
 * For each count in number of particles to release
 * -# Keep ready a particle with randomized attributes.
 * -# If lostParticle bound is not 0, then we still have more indices to overwrite, so insert it into that.
 * -# If not, then simply insert particles at the end of the array.
 * 
 * \note End of the array does NOT correspond to end of the domain. There is no relation between the array position
 * and position inside the domain.
 *
 */
void releaseParticles(
        const int numParticlesToRelease,  /**< Number of particles to release */
        const Particle* inputData,        /**< Input data of MD simulated particles to choose from */
        const int inputCount,             /**< Size of the MD input data */
        Particle* domainParticles,        /**< An array with the Particles currently in the domain */
        int* domainParticleCount,         /**< Count of particles in the domain */
        int* lostParticlesArray,          /**< Array indicating the particle indices where the particles have left the domain */
        int* lostParticleBound,            /**< Pass by pointer for modifying lostParticles index counts */
        unsigned int *randSeeds
        )
{
    int tid = omp_get_thread_num();

    int i;

    // TODO: move this one level above?
    // increment the domainParticle count in preparation for insertion
    //(*domainParticleBound) += (numParticlesToRelease - (*lostParticlesBound));
    //int particleBoundCopy = (*domainParticleBound);

    // assert just in case
    //assert(*domainParticleCount <= PARTICLE_SIZE);

    // loop through MD_data and randomize particle attributes
    #pragma omp for schedule(static)
    for(i = 0; i < numParticlesToRelease; i++)
    {
        Particle releasedParticle = randomizeParticleAttribs(inputData[rand_r( &(randSeeds[tid]) ) % inputCount] );
        //Particle releasedParticle = randomizeParticleAttribs(inputData[rand() % inputCount] );

        // if bound is not 0
        // TODO: avoid this somehow?!
        #pragma omp critical
        {
        if(*lostParticleBound >= 0)
        {
            // insert particles from the end
            domainParticles[ lostParticlesArray[(*lostParticleBound)] ] = releasedParticle;
            (*lostParticleBound)--;
            //(*domainParticleBound)--;
        }
        // if it is 0
        // then insert remaining elements at the end of the domainParticle array
        else
        {
            // store a copy of domainParticleBound
            domainParticles[ *domainParticleCount] = releasedParticle;

            // increment counter since we are inserting at the end of the
            // domainParticles array
            (*domainParticleCount)++;
        }
        } // end of pragma critical
    }

}

/*!
 * If when this function is called, lostParticlesBound is still not zero,
 * that means we lost a large enough number of particles, that the
 * releasedParticles did not fill up the gaps.
 * So swap all gaps with the last few particles and decrement bounds
 * appropriately
 */
void swapGapsWithEndParticles(Particle* domainParticles, int* domainParticleCount,
                              int* lostParticlesArray, int* lostParticleBound)
{
    int i, j = 0;
    int endIndex = (*domainParticleCount)-1;

    // TODO: looks like it will be hard to parallleize this
    // since threads can see the same domainParticleCount and end up
    // copying the same particle into two different places!
    #pragma omp for schedule(static)
    for(i = (*lostParticleBound); i >= 0; i--)
    {
        j = (*lostParticleBound)-i;
        // swap lostParticlesBound with domainParticleBound value
        domainParticles[ lostParticlesArray[i] ] = domainParticles[ endIndex-j ];
    }

    #pragma omp single
    {
    (*domainParticleCount) -= (*lostParticleBound)+1;   // loop executes (*lostPBound)+1 times
    (*lostParticleBound) = -1;
    }


    // at the end of this, lostParticleBound HAS to be -1 again no?
    //assert(*lostParticleBound == -1);
}

int moveParticlesInField(Particle* domainParticles, int domainParticleBound,
                         int* lostParticlesArray,//int* lostParticlesBound,
                         EField* ElectricField, GridInfo* gInfo)
{
    const double halfTimeStep = 0.5 * T_PIC;

    // declare arrays which hold the weighted values
    double weights[8] = {0};
    //EField  efield[8];          // how to zero initialize this by default?
    int    indices[8][3];
    EField elecField;

    // loop through all particles
    int i, lostParticleCount = 0;
    #pragma omp for schedule(static)
    for(i = 0; i < domainParticleBound; i++)
    {
        Particle* p = &domainParticles[i];

        // TODO: optimize this by avoiding reweighting corner efield values
        // if the new particle is in the same cell as the previous?
        // get the cellWeights
        getCellWeights(p->x, p->y, p->z,
                       gInfo->spacing, gInfo->invSpacing,
                       indices, weights);

        // reset the elecfield components to zero
        // in preparation for weighting
        elecField.components[0] = 0;
        elecField.components[1] = 0;
        elecField.components[2] = 0;

        // now store the ElectricField value
        int c, j;
        for(c = 0; c < 8; c++)
        {
            int index = INDEX_1D(gInfo->numNodes, indices[c][0], indices[c][1], indices[c][2]);
            double weightFactor = weights[c];
            EField efield = ElectricField[ index ];
            for(j = 0; j < 3; j++)
                elecField.components[j] += weightFactor * efield.components[j];
        }

        // TODO: generalize as a loop through the number of components?
        const double multFactor = (p->charge * ELECTRONIC_CHARGE * T_PIC/ p->mass);

        // update the velocity
        // save the old velocities
        double old_Vx = p->Vx,
               old_Vy = p->Vy,
               old_Vz = p->Vz;

        p->Vx += elecField.components[0] * multFactor;
        p->Vy += elecField.components[1] * multFactor;
        p->Vz += elecField.components[2] * multFactor;

        // update the position
        p->x += halfTimeStep * (old_Vx + p->Vx);
        p->y += halfTimeStep * (old_Vy + p->Vy);
        p->z += halfTimeStep * (old_Vz + p->Vz);

        // check if the particle is out of bounds
        if(!isParticleInDomain(*p))
        {
            // if so, increment lost particle counter
            lostParticlesArray[lostParticleCount] = i;
            lostParticleCount++;
        }
    }

    return lostParticleCount;
}

int simulateSingleParticleTillDomainExit(
        const Particle *input_data,
        const int inputCount,
        EField *ElectricField,
        GridInfo *gridInfo
        )
{
    Particle p;
    int lostParticles[1];
    int particleCount = 1;
    int lostParticleBound = -1;

    // set the pArr to the particle with the least X-Velocity
    // from the input set
    p.Vx = 1e10;
    int i;
    for(i = 0; i < inputCount; i++)
    {
        if((input_data[i].Vx < p.Vx) && (input_data[i].Vx >= 0) )
            p = input_data[i];
    }

    // release a single particle
    //releaseParticles(1,
    //                 input_data, inputCount,
    //                 pArr, &particleCount,
    //                 lostParticles,
    //                 &lostParticleBound);

    // count the number of iterations required till steady state
    int possibIterCount = 0;
    while(isParticleInDomain(p))
    {
        moveParticlesInField(&p, particleCount,
                             lostParticles,//&lostParticleBound,
                             ElectricField, gridInfo);
        possibIterCount++;
    }

    return possibIterCount;
}


// routine to re-sort particles based on their distance
// from capillary center
void resortParticles(Particle *domainParticles, const int particleCount)
{
    qsort(domainParticles, particleCount, sizeof(Particle), axial_dist_compare);
}

#endif
