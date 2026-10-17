# Case A1 – Hyperbolic PDE Upwinding

**File:** `A1_hyperbolic_pde_upwinding.m`  
---

## 1. The Physics

We solve the simplest hyperbolic conservation law – the **linear advection equation**:

```
∂u/∂t + a * ∂u/∂x = 0 ,   a > 0
```

- Domain: periodic interval [0, 1]
- Wave speed: a = 1
- Initial condition: a square wave (discontinuity) + a smooth Gaussian pulse

The exact solution is simply the initial profile shifted by the distance a*t:

```
u_exact(x,t) = u0(x - a*t)
```

(with wrap-around because of periodicity).

---

## 2. Why we stop at t = 1.0

```
a * t = 1 * 1 = 1
```

The wave travels **exactly one full domain length**.  
Therefore the exact solution at t = 1 is identical to the initial condition (just shifted by a whole period).

This gives a perfect visual test:
- Black dashed line = initial condition
- Black solid line = exact solution at t = 1 (should lie on top of the initial condition)
- Any deviation of the coloured curves from the black line is pure numerical error

---

## 3. Key parameters

| Parameter | Value | Meaning |
|-----------|-------|---------|
| Nx        | 200   | Number of cells (grid resolution).  dx = 1/200 = 0.005 |
| CFL       | 0.4   | Courant–Friedrichs–Lewy number (see below) |
| a         | 1.0   | Advection speed |
| Tfinal    | 1.0   | Final time = one full period |

### What is the CFL number?

```
CFL = (a * dt) / dx
```

It is a dimensionless stability limit for explicit schemes.  
Physically it means: “In one time step the wave is allowed to travel at most a fraction CFL of one cell.”

- Most explicit upwind schemes require **CFL ≤ 1** for stability.
- We choose 0.4 as a comfortable safety margin.
- The time step is then computed as:  dt = CFL * dx / a

---

## 4. The three numerical schemes

### 4.1 First-order Upwind

**Idea**  
Information travels only in the direction of the “wind”.  
Because a > 0, the value that enters cell i comes from the left neighbour i-1.

**Update formula**

```
u_i^{n+1} = u_i^n - (a*dt/dx) * (u_i^n - u_{i-1}^n)
```

**Properties**
- Extremely robust and monotone (never creates new maxima or minima)
- Only first-order accurate → highly dissipative
- The square wave is turned into a smooth ramp after one period

This is the simplest possible upwind scheme and the one most CFD codes fall back to when everything else fails.

### 4.2 Lax–Friedrichs

**Idea**  
Replace the solution at the new time level by an average of the two neighbours, then add a central difference for the flux.  
The averaging introduces a large amount of numerical viscosity.

**Update formula** (classic form)

```
u_i^{n+1} = (u_{i-1}^n + u_{i+1}^n)/2  -  (a*dt)/(2*dx) * (u_{i+1}^n - u_{i-1}^n)
```

**Properties**
- Also first-order and monotone under CFL ≤ 1
- Even more dissipative than plain upwind
- Very easy to code, historically important, but rarely used in production today

### 4.3 MUSCL + minmod limiter

**MUSCL** = Monotone Upstream-centered Scheme for Conservation Laws (van Leer, 1979)

**Basic idea**  
Instead of assuming the solution is constant inside each cell (first-order), we reconstruct a **linear** profile:

```
u(x) ≈ u_i + σ_i * (x - x_i)
```

We then evaluate this profile at the cell faces to obtain a more accurate flux.

**The problem**  
A pure unlimited slope produces oscillations near discontinuities (Godunov’s theorem).  
We therefore **limit** the slope.

**minmod limiter**

```
σ_i = minmod( u_i - u_{i-1} ,  u_{i+1} - u_i )
```

where

```
minmod(a,b) =
    a     if |a| < |b|  and  a*b > 0
    b     if |b| < |a|  and  a*b > 0
    0     otherwise
```

- Opposite signs → local extremum → set slope to zero (locally first-order)
- Same sign → take the smaller difference → the scheme stays monotone (TVD)

**Result**
- Second-order accurate in smooth regions (the Gaussian stays almost perfect)
- Still monotone at discontinuities (the square wave remains sharp, only the corners are slightly rounded)

This combination is why limited MUSCL (and its modern relatives WENO, PPM, etc.) is the everyday work-horse of industrial CFD codes.

---

## 5. Expected visual result

| Scheme              | Order        | Dissipation | Oscillations | Appearance after one period             |
|---------------------|--------------|-------------|--------------|-----------------------------------------|
| Exact (black)       | —            | None        | None         | Perfect square + Gaussian               |
| 1st-order Upwind    | 1            | High        | None         | Square becomes a smooth ramp            |
| Lax–Friedrichs      | 1            | Very high   | None         | Even more smeared                       |
| MUSCL + minmod      | 2 (smooth)   | Low         | None         | Sharp square, almost exact Gaussian     |

The green MUSCL curve should lie **between** the exact solution and the first-order upwind result: sharper than upwind, no oscillations.

---

