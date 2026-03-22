# Physics Simulations

A collection of classical physics simulations built in Python using NumPy and Matplotlib. Started as a transition from MATLAB, these scripts explore numerical methods and animation.

---

## Simulations

### Projectile Motion (No Drag)
`projectile_no_drag.py`

Analytical solution for idealized projectile motion. Computes the full trajectory in one shot using vectorized NumPy operations — no loop needed.

- Analytical solution (vectorized)
- Configurable launch speed and angle
- Animated ball along trajectory

### Projectile Motion (With Drag)
`projectile_drag.py`

Adds quadratic air resistance, which makes an analytical solution impossible. Uses the Euler method for numerical integration.

- Quadratic drag model: `F = -b|v|v`
- Euler method numerical integration
- Ground condition with automatic time window
- Animated ball along trajectory

### Real Pendulum
`real_pendulum.py`

Driven damped pendulum simulation using the nonlinear equation of motion. Uses the Euler-Cromer method to conserve energy.

- Nonlinear `sin(θ)` equation of motion
- Damping and external driving force
- Euler-Cromer integration
- 3-panel live animation: pendulum, θ/ω vs time, phase space

---

## Numerical Methods Used

| Method | Used In | Notes |
|--------|---------|-------|
| Vectorized (analytical) | Projectile no drag | No loop needed |
| Euler | Projectile with drag | Simple but energy drift |
| Euler-Cromer | Real pendulum | Better energy conservation for oscillators |

---

## Requirements

```
pip install numpy matplotlib
```

## Running

```bash
python projectile_no_drag.py
python projectile_drag.py
python real_pendulum.py
```

---

## Roadmap

- [ ] Wave equation (1D string)
- [ ] 2D drum membrane
- [ ] Monte Carlo simulation
- [ ] N-body gravity

# Gyroscope Animation

This program visualizes the phenomenon of gyroscopic behavior using vector dynamics.

- **Linear vectors:** $\vec{r}$, $\vec{p}$, $\vec{F}$
- **Angular vectors:** $\vec{r}$, $\vec{L}$, $\vec{\tau}$

The simulation demonstrates a spinning object with large angular momentum. In the absence of external dissipation (e.g., friction), the system exhibits ideal gyroscopic precession, where the angular momentum vector changes direction while maintaining its magnitude.

Friction and other damping effects are neglected to allow the gyroscopic motion to persist indefinitely and clearly illustrate the underlying physics.

### Vector Definitions

- $\vec{r}$: Position vector from the axis of rotation to the point of application of the force $\vec{F}$

- $\vec{p}$: Linear momentum vector, defined as  
  $\vec{p} = m\vec{v}$, and pointing in the same direction as the velocity $\vec{v}$

- $\vec{F}$: Force vector, defined by  
  $\vec{F} = \frac{d\vec{p}}{dt}$

- $\vec{L}$: Angular momentum vector, defined as  
  $\vec{L} = \vec{r} \times \vec{p}$  
  It is perpendicular to the plane formed by $\vec{r}$ and $\vec{p}$

- $\vec{\tau}$: Torque vector, defined as  
  $\vec{\tau} = \vec{r} \times \vec{F}$  
  It is perpendicular to the plane formed by $\vec{r}$ and $\vec{F}$, and satisfies  
  $\vec{\tau} = \frac{d\vec{L}}{dt}$