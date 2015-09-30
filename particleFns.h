#ifndef PARTICLEFNS_H
#define PARTICLEFNS_H


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
    //Particle ret;                                                                       

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

/*!
 * This routine releases the a specified number of particles into the domain. It uses the lostParticles index
 * count array to keep track of indices of particles which have left the domain. The outline is as follows:
 * For each count in number of particles to release
 * -# Keep ready a particle with randomized attributes.
 * -# If lostParticle bound is not 0, then we still have more indices to overwrite, so insert it into that.
 * -# If not, then simply insert particiles at the end of the array.
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


#endif
