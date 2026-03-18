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
