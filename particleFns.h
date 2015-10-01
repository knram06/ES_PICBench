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

double _getSquaredDistance(double x1, double y1, double z1, double x2, double y2, double z2)
{
    double tx = (x1 - x2);
    double ty = (y1 - y2);
    double tz = (z1 - z2);

    return (tx*tx + ty*ty + tz*tz);
}

int updateChargeFractions(
        const Particle *particleList, int particleCount,
        double *chargeFractionList,
        int *chargeFractionIndices,
        GridInfo *gInfo
        )
{
    int i, validNodesNumber = 0;
    int searchIndex = 0;

    int numNodes = gInfo->numNodes;
    double invSpacing = gInfo->invSpacing;
    double spacing = gInfo->spacing;

    // cache the multiplication factor
    // taking spacing factors to rhs and taking vol as spacing^3
    const double multFactor = -invSpacing * ELECTRONIC_CHARGE / FREE_SPACE_PERMITTIVITY;

    for(i = 0; i < particleCount; i++)
    {
        const Particle *p = &particleList[i];

        // round off and get the nearest indices - these must be the lower
        // corner indices of a cell - due to nature of rounding?
        int iPos = (int)(round( (p->x)*invSpacing) );
        int jPos = (int)(round( (p->y)*invSpacing) );
        int kPos = (int)(round( (p->z)*invSpacing) );

        // array to store calculated distances from vertices
        // distances can be squared as it does not affect the comparison
        double distances[8];
        int indices[8][3];

        // manually fill in the distances and indices to each of the vertices
        int ti = iPos, tj = jPos, tk = kPos;
        indices[0][0] = ti; indices[0][1] = tj; indices[0][2] = tk;
        distances[0] = _getSquaredDistance(p->x, p->y, p->z, ti*spacing, tj*spacing, tk*spacing);

        // even if it is redundant, explicitly setting to make it clear
        ti = iPos; tj = jPos; tk = kPos+1;
        indices[1][0] = ti; indices[1][1] = tj; indices[1][2] = tk;
        distances[1] = _getSquaredDistance(p->x, p->y, p->z, ti*spacing, tj*spacing, tk*spacing);

        ti = iPos; tj = jPos+1; tk = kPos;
        indices[2][0] = ti; indices[2][1] = tj; indices[2][2] = tk;
        distances[2] = _getSquaredDistance(p->x, p->y, p->z, ti*spacing, tj*spacing, tk*spacing);

        ti = iPos; tj = jPos+1; tk = kPos+1;
        indices[3][0] = ti; indices[3][1] = tj; indices[3][2] = tk;
        distances[3] = _getSquaredDistance(p->x, p->y, p->z, ti*spacing, tj*spacing, tk*spacing);

        /*********************************************/
        /*********************************************/
        ti = iPos+1, tj = jPos, tk = kPos;
        indices[4][0] = ti; indices[4][1] = tj; indices[4][2] = tk;
        distances[4] = _getSquaredDistance(p->x, p->y, p->z, ti*spacing, tj*spacing, tk*spacing);

        // even if it is redundant, explicitly setting to make it clear
        ti = iPos+1; tj = jPos; tk = kPos+1;
        indices[5][0] = ti; indices[5][1] = tj; indices[5][2] = tk;
        distances[5] = _getSquaredDistance(p->x, p->y, p->z, ti*spacing, tj*spacing, tk*spacing);

        ti = iPos+1; tj = jPos+1; tk = kPos;
        indices[6][0] = ti; indices[6][1] = tj; indices[6][2] = tk;
        distances[6] = _getSquaredDistance(p->x, p->y, p->z, ti*spacing, tj*spacing, tk*spacing);

        ti = iPos+1; tj = jPos+1; tk = kPos+1;
        indices[7][0] = ti; indices[7][1] = tj; indices[7][2] = tk;
        distances[7] = _getSquaredDistance(p->x, p->y, p->z, ti*spacing, tj*spacing, tk*spacing);

        // MIN NODE DISTANCE
        // FOR CHARGE ASSIGNMENT - for now
        // now get the minindex corresponding to minDistance
        int j;
        double minDistance = 1e8;
        int minIndex;
        for(j = 0; j < 8; j++)
        {
            if(distances[j] < minDistance)
            {
                minDistance = distances[j];
                minIndex = j;
            }
        } // end of for loop for min distance

        // now if the point is in the interior
        // then store it
        // Why do we do this? Since points on the boundaries have the BCs imposed,
        // we don't want to overwrite it in the Solver data
        // so we specify an array of indices explicitly
        ti = indices[minIndex][0];
        tj = indices[minIndex][1];
        tk = indices[minIndex][2];

        // if within interior points of the grid
        if( ((ti > 0) && (ti < numNodes-1))
                      &&
            ((tj > 0) && (tj < numNodes-1))
                      &&
            ((tk > 0) && (tk < numNodes-1)) )
        {
            int pos = INDEX_1D(numNodes, ti, tj, tk);

            // wrap around index
            // reusing last searched value to speed things up?
            for(j = 0; j < validNodesNumber; j++)
            {
                searchIndex = searchIndex % validNodesNumber;
                if(chargeFractionIndices[searchIndex] == pos)
                    chargeFractionList[searchIndex] += p->charge * multFactor;
            }

            // if pos was not found, then write into and expand size
            if(j == validNodesNumber)
            {
                chargeFractionIndices[validNodesNumber] = pos;
                chargeFractionList[validNodesNumber] = p->charge * multFactor;
                validNodesNumber++;
            }
        } // end of if check for indices within valid limits
        
    } // end of loop through particles array

    return validNodesNumber;
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
        int* lostParticleBound            /**< Pass by pointer for modifying lostParticles index counts */
        )
{
    int i;

    // TODO: move this one level above?
    // increment the domainParticle count in preparation for insertion
    //(*domainParticleBound) += (numParticlesToRelease - (*lostParticlesBound));
    //int particleBoundCopy = (*domainParticleBound);

    // assert just in case
    assert(*domainParticleCount <= PARTICLE_SIZE);

    // loop through MD_data and randomize particle attributes
    for(i = 0; i < numParticlesToRelease; i++)
    {
        Particle releasedParticle = randomizeParticleAttribs(inputData[rand() % inputCount] );

        // if bound is not 0
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
    while(*lostParticleBound >= 0)
    {
        // swap lostParticlesBound with domainParticleBound value
        domainParticles[ lostParticlesArray[*lostParticleBound] ] = domainParticles[ (*domainParticleCount) - 1 ];
        (*domainParticleCount)--;
        (*lostParticleBound)--;
    }

    // at the end of this, lostParticleBound HAS to be -1 again no?
    //assert(*lostParticleBound == -1);
}

int moveParticlesInField(Particle* domainParticles, int domainParticleBound,
                         int* lostParticlesArray,//int* lostParticlesBound,
                         EField* ElectricField, GridInfo* gInfo)
{
    const double invSpacing = gInfo->invSpacing;
    const double halfTimeStep = 0.5 * T_PIC;

    // FOR NOW
    // Approximate each particle to nearest node
    // TODO: Weight it based on distances to other nodes?

    // loop through all particles
    int i, lostParticleCount = 0;
    //int lostParticlesBound = -1;
    for(i = 0; i < domainParticleBound; i++)
    {
        Particle* p = &domainParticles[i];

        // round off and get the nearest indices
        int iPos = (int)(round( (p->x)*invSpacing) );
        int jPos = (int)(round( (p->y)*invSpacing) );
        int kPos = (int)(round( (p->z)*invSpacing) );

        // now store the ElectricField value
        EField* elecField = &GRID_1D(ElectricField, iPos, jPos, kPos);

        // TODO: generalize as a loop through the number of components?
        const double multFactor = (p->charge * ELECTRONIC_CHARGE * T_PIC/ p->mass);

        // update the velocity
        // save the old velocities
        double old_Vx = p->Vx,
               old_Vy = p->Vy,
               old_Vz = p->Vz;

        p->Vx += elecField->components[0] * multFactor;
        p->Vy += elecField->components[1] * multFactor;
        p->Vz += elecField->components[2] * multFactor;

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

// routine to re-sort particles based on their distance
// from capillary center
void resortParticles(Particle *domainParticles, const int particleCount)
{
    qsort(domainParticles, particleCount, sizeof(Particle), axial_dist_compare);
}

#endif
