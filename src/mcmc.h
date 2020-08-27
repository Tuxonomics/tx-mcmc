// mcmc.h
//
// A lot of the following code was influenced by Stan and can thus be seen as
// a derivative work. Stan is subject to the BSD 3-clause license.
//
// BSD 3-Clause License
//
// Copyright (c) 2011-2020, Stan Developers and their Assignees
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its contributors
//   may be used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
// OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
// IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include "utils.h"
#include "../deps/prng.h"
#include "matrix.h"
#include "mvnormal.h"


#define AL _ALLOCATOR_DEFAULT


// -------------------- //
// State of Hamiltonian //
// -------------------- //

struct HamiltonianState {
    f64    energy;
    Matrix q;
    Matrix p;
};


void HamiltonianStateInit(
    Allocator al,
    HamiltonianState *hs,
    u32 n
) {
    hs->energy = -INFINITY;

    MInit( al, &(hs->q), n, 1 );
    MInit( al, &(hs->p), n, 1 );
}


void HamiltonianStateFree(
    Allocator al,
    HamiltonianState *hs
) {
    hs->energy = NAN;
    MFree( al, &(hs->q) );
    MFree( al, &(hs->p) );
}


void HamiltonianStateCopy(
    HamiltonianState *lhs,
    HamiltonianState rhs
) {
    ASSERT_MSG( lhs->q.rows == rhs.q.rows && lhs->q.cols == rhs.q.cols,
                "HamiltonianStateCopy: q in lhs and rhs need to be of equal size." );
    ASSERT_MSG( lhs->p.rows == rhs.p.rows && lhs->p.cols == rhs.p.cols,
                "HamiltonianStateCopy: q in lhs and rhs need to be of equal size." );
    lhs->energy = rhs.energy;
    MCopy( lhs->q, rhs.q );
    MCopy( lhs->p, rhs.p );
}



// ----------- //
// Hamiltonian //
// ----------- //

template <typename F>
INLINE
f64 PotentialEnergy( HamiltonianState hs, F f ) {
    return f( hs.q.data, hs.q.rows );
}


f64 KineticEnergy( HamiltonianState hs, Matrix iMass )
{
    ASSERT_MSG( hs.q.rows == hs.p.rows && hs.q.cols == 1 && hs.p.cols == 1,
                "KineticEnergy: q and p need to be of equal size." );

    u32 n = hs.q.rows;

    f64 res = 0.0;

    Mat(tmp0);
    Mat(tmp1);

    MInit( AL, &tmp0, n, 1 );
    MInit( AL, &tmp1, 1, 1 );
    defer({
        MFree( AL, &tmp0 );
        MFree( AL, &tmp1 );
    });

    MMul( tmp0, iMass, hs.p );

    Mat(pT);
    pT.data = hs.p.data; pT.rows = 1; pT.cols = n;

    MMul( tmp1, pT, tmp0 );

    return -0.5 * tmp1[0];
}


template <typename F>
INLINE
f64 Hamiltonian( HamiltonianState hs, Matrix iMass, F f ) {
    return PotentialEnergy( hs, f ) + KineticEnergy( hs, iMass );
}


f64 acceptanceRatio( HamiltonianState os, HamiltonianState ns )
{
    if ( ! isnan(ns.energy) ) {
        return exp( ns.energy - os.energy );
    }
    return 0.0;
}


// ---------- //
// Integrator //
// ---------- //

template <typename G>
void LeapfrogStep(
    HamiltonianState os, HamiltonianState ns, Matrix iMass, f64 eps, i32 dir, G gradient
) {
    ASSERT_MSG( os.q.rows == iMass.rows && os.q.rows == iMass.cols,
                "LeapfrogStep: dimensions of q and iMass don't match." );

    f64 step = dir * eps;

    u32 n = os.q.rows;

    Mat(tmp0); Mat(tmp1); Mat(tmp2);
    MInit(AL, &tmp0, n, 1);
    MInit(AL, &tmp1, n, 1);
    MInit(AL, &tmp2, n, 1);
    defer({
        MFree(AL, &tmp0);
        MFree(AL, &tmp1);
        MFree(AL, &tmp2);
    });

    Mat(grad);
    MInit(AL, &grad, n, 1);
    defer({
        MFree(AL, &grad);
    });

    gradient(grad.data, os.q.data, os.q.rows);

    MScale(tmp0, grad, 0.5*step);
    MAdd(tmp1, tmp0, os.p);

    MMul(tmp0, iMass, tmp1);
    MScale(tmp2, tmp1, step);
    MAdd(ns.q, tmp1, tmp2);

    gradient(grad.data, ns.q.data, ns.q.rows);

    MScale(tmp2, grad, 0.5*step);
    MAdd(ns.p, tmp2, tmp1);
}


// -------------- //
// MCMC Utilities //
// -------------- //

#define N_DIAGNOSTICS 5

// NOTE(jonas): diagnostics save
// - acceptance ratio
// - if chain is divergent
// - tree depth for NUTS
// - number of steps taken
// - energy of Hamiltonian

struct Results {
    Matrix sample;
    Matrix diagnostics;
};


void ResultsInit( Allocator al, Results *r, u32 nIter, u32 nParams )
{
    MInit( al, &(r->sample), nIter, nParams );
    MInit( al, &(r->diagnostics), nIter, N_DIAGNOSTICS);
}


void ResultsFree( Allocator al, Results *r )
{
    MFree( al, &(r->sample) );
    MFree( al, &(r->diagnostics) );
}


#undef N_DIAGNOSTICS




// ----------------------- //
// Hamiltonian Monte Carlo //
// ----------------------- //

template <typename F, typename G>
void HMCStep(
    Results r,
    prng::Generator rng,
    f64 eps, Matrix iMass, Matrix iMassCT, u32 it, u32 nSteps,
    F f, G gradient
) {
    Mat(q0);
    q0.data = r.sample.data + (it-1) * r.sample.cols;
    q0.rows = r.sample.cols; q0.cols = 1;

    u32 n = q0.rows;

    Mat(p0);
    MInit( AL, &p0, n, 1 );
    defer({ MFree(AL, &p0); });

    HamiltonianState os;
    HamiltonianStateInit( AL, &os, n );
    HamiltonianState ns;
    HamiltonianStateInit( AL, &ns, n );

    defer({
        HamiltonianStateFree( AL, &os );
        HamiltonianStateFree( AL, &ns );
    });

    ASSERT_MSG( ns.q.rows == q0.rows, "HMCStep: ns.q and q0 don't match." );
    MCopy(ns.q, q0);
    MvNormal(rng, iMassCT, p0);

    ASSERT_MSG( ns.p.rows == p0.rows, "HMCStep: ns.p and p0 don't match." );
    MCopy(ns.p, p0);

    for ( u32 i = 0; i < nSteps; ++i ) {
        HamiltonianStateCopy( &os, ns );
        LeapfrogStep(os, ns, iMass, eps, 1, gradient);
    }

    ns.energy = Hamiltonian(ns, iMass, f);

    HamiltonianState tmpState = { 0 };
    tmpState.q = q0; tmpState.p = p0;
    tmpState.energy = Hamiltonian( tmpState, iMass, f );

    f64 a = acceptanceRatio(tmpState, ns);
    f64 u = prng::Uniform(rng);

    if ( (a < 1.0) && (u > a) ) {
        MSetRow(r.sample, it, q0);
        // add diagnostics
    }
    else {
        MSetRow(r.sample, it, ns.q);
        // add diagnostics
    }
}


template <typename F, typename G>
Results RunHMC(
    Allocator al,
    f64       eps,
    u32       intTime,
    Matrix    iMass,
    u32       nIter,
    Matrix    initialValues,
    u64       seed,
    F         f,
    G         gradient
) {
    ASSERT( initialValues.rows == 1 || initialValues.cols == 1 );

    prng::Xoshiro256StarStar state;
    prng::Generator rng = prng::InitXoshiro256StarStar( &state, seed );

    u32 nParams = initialValues.rows*initialValues.cols;

    u32 nSteps = u32((f64)intTime / eps);
    nSteps = ( nSteps < 1 ) ? 1 : nSteps;

    // Prepare inverse mass matrix for normal draws
    Mat(iMassCT);
    MInit(AL, &iMassCT, iMass.rows, iMass.cols);
    Mat(tmp);
    MInit(AL, &tmp, iMass.rows, iMass.cols);
    defer({
        MFree(AL, &iMassCT);
        MFree(AL, &tmp);
    });

    MCopy(tmp, iMass);
    MCholeski(tmp, 'U');
    MTranspose(iMassCT, tmp);


    Results r = { 0 };
    ResultsInit( al, &r, nIter, nParams );

    MSetRow( r.sample, 0, initialValues );

    for ( u32 it = 1; it < nIter; ++it ) {
        HMCStep(r, rng, eps, iMass, iMassCT, it, nSteps, f, gradient);
    }

    return r;
}


#undef AL

