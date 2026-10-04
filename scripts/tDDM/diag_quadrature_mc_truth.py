# Monte Carlo ground truth for the two-phase diffusion (unit diffusion, sigma-rescaled
# units as in the engine): start bsep/2, drift d1 before tau and d2 after, absorbing at
# 0 (choice_ddm = 0) and bsep. Euler steps of 5e-5 s with a Brownian-bridge crossing
# correction. Density of a lower-boundary hit estimated in windows around dt.
import numpy as np, pandas as pd, time
cases = pd.read_csv("mc_cases.csv")
rng = np.random.default_rng(20260927)
N, H = 400_000, 5e-5
for _, c in cases.iterrows():
    t0 = time.time()
    T = c["dt"] + 0.08
    x = np.full(N, c.bsep / 2); alive = np.ones(N, bool); hit = np.full(N, np.nan)
    t = 0.0; sq = np.sqrt(H)
    while t < T and alive.any():
        v = c.d1 if t < c.tau else c.d2
        idx = np.flatnonzero(alive)
        xo = x[idx]; xn = xo + v * H + sq * rng.standard_normal(idx.size)
        # bridge crossing probability for the lower boundary between xo and xn
        pl = np.where((xo > 0) & (xn > 0), np.exp(-2 * xo * xn / H), 1.0)
        pu = np.where((xo < c.bsep) & (xn < c.bsep), np.exp(-2 * (c.bsep - xo) * (c.bsep - xn) / H), 1.0)
        u = rng.random(idx.size)
        low = u < pl; up = (~low) & (rng.random(idx.size) < pu)
        hit[idx[low]] = t + H / 2
        alive[idx[low | up]] = False
        x[idx] = xn
        t += H
    def dens(a, b): return np.mean((hit >= a) & (hit < b)) / (b - a)
    print(f"subject {int(c.subject)}  tau {c.tau:.4f}  dt {c["dt"]:.4f}  ({time.time()-t0:.0f}s)")
    print(f"  MC density  [dt-20ms, dt)     {dens(c["dt"]-0.02, c["dt"]):.4f}")
    print(f"  MC density  [dt, dt+20ms)     {dens(c["dt"], c["dt"]+0.02):.4f}")
    print(f"  MC density  [dt+40, dt+60ms)  {dens(c["dt"]+0.04, c["dt"]+0.06):.4f}")
    print(f"  engine grid30 at dt           {c.grid30:.4f}")
    print(f"  engine grid30 at dt+50ms      {c.grid30_plus50ms:.4f}")
