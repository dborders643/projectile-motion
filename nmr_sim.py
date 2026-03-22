"""
NMR Hydrogen Spin Ensemble Simulation
======================================
Shows N hydrogen magnetization vectors on a Bloch sphere, precessing
and dephasing in real time after an RF pulse.

Physics (rotating frame, scaled time):
  - Each spin i has offset frequency delta_i drawn from a Lorentzian
    distribution (field inhomogeneity / T2* effect)
  - After a flip-angle RF pulse, all spins tip into the XY plane
  - They precess at delta_i and relax: Mxy *= exp(-t/T2), Mz -> 1
  - Net FID = sum of all Mxy vectors → FFT gives the spectrum

Panels:
  [Bloch sphere 3D]   [Top-down XY precession]
  [FID time signal]   [FFT frequency spectrum]

Controls (left panel):
  T1, T2, flip angle, inhomogeneity (ppm), N spins, speed
  Buttons: RF Pulse, Reset
  SPACE = apply RF pulse
  Right-click drag on Bloch sphere = rotate view

Dependencies:
  pip install numpy pygame scipy

Run:
  python nmr_simulation.py
"""

import math
import numpy as np
import pygame
from scipy.fft import fft, fftfreq

# ── Constants ────────────────────────────────────────────────────────────────
GAMMA  = 267.5e6        # proton gyromagnetic ratio rad/T/s
TWO_PI = 2 * math.pi
DT     = 5e-4           # simulation time step (s)
T_MAX  = 0.6            # total signal accumulation time (s)

# ── Colours ───────────────────────────────────────────────────────────────────
BG          = ( 12,  12,  30)
PANEL_BG    = ( 18,  18,  46)
CTRL_BG     = ( 22,  22,  54)
ACCENT      = (  0, 212, 170)
ACCENT2     = (255, 180,  60)
WHITE       = (220, 220, 240)
GREY        = (100, 100, 140)
DIM         = ( 40,  40,  80)
RED_AX      = (255,  80,  80)
GREEN_AX    = ( 80, 220, 100)
YELLOW_AX   = (255, 210,  60)
SLIDER_TRK  = ( 35,  35,  75)
SLIDER_FILL = (  0, 180, 140)
BTN_NRM     = ( 38,  38,  88)
BTN_HOV     = ( 58,  58, 128)
BTN_ACT     = (  0, 160, 120)


# ── Physics: Spin Ensemble ────────────────────────────────────────────────────
class SpinEnsemble:
    """
    N hydrogen spins in the rotating frame.

    Each spin has a slightly different precession frequency (delta_omega)
    drawn from a Lorentzian distribution, modelling B0 field inhomogeneity.
    After an RF pulse tips the spins, they fan out in the XY plane as they
    precess at their individual offset frequencies — this dephasing is what
    causes the FID envelope to decay (T2* effect).
    """

    def __init__(self, N=24, T1=1.0, T2=0.08, dB0_ppm=5.0, flip_deg=90.0):
        self.N        = N
        self.T1       = T1
        self.T2       = T2
        self.dB0_ppm  = dB0_ppm
        self.flip_deg = flip_deg
        self.t        = 0.0
        self._rng     = np.random.default_rng(42)
        self._build_offsets()
        self._init_equilibrium()

    def _build_offsets(self):
        # Lorentzian (Cauchy) frequency spread in Hz
        sigma = self.dB0_ppm * GAMMA * 3.0 / 1e6 / TWO_PI
        self.delta_omega = self._rng.standard_cauchy(self.N) * sigma  # rad/s
        # Clamp outliers so display stays sensible
        limit = sigma * 12
        self.delta_omega = np.clip(self.delta_omega, -limit, limit)

    def _init_equilibrium(self):
        self.Mx = np.zeros(self.N)
        self.My = np.zeros(self.N)
        self.Mz = np.ones(self.N)
        self.t  = 0.0
        self._fid_t   = []
        self._fid_sig = []

    def apply_rf_pulse(self):
        """Tip all spins by flip_deg around the x-axis."""
        flip = math.radians(self.flip_deg)
        cf, sf = math.cos(flip), math.sin(flip)
        new_My =  self.My * cf + self.Mz * sf
        new_Mz = -self.My * sf + self.Mz * cf
        self.My = new_My
        self.Mz = new_Mz

    def step(self, dt):
        """Advance simulation by dt seconds."""
        # Precession at each spin's offset frequency
        phase = self.delta_omega * dt
        cp, sp = np.cos(phase), np.sin(phase)
        new_Mx = self.Mx * cp - self.My * sp
        new_My = self.Mx * sp + self.My * cp
        self.Mx, self.My = new_Mx, new_My

        # T2 transverse relaxation
        d2 = math.exp(-dt / self.T2)
        self.Mx *= d2
        self.My *= d2

        # T1 longitudinal recovery
        d1 = math.exp(-dt / self.T1)
        self.Mz = 1.0 - (1.0 - self.Mz) * d1

        self.t += dt
        self._fid_t.append(self.t)
        self._fid_sig.append(complex(float(np.mean(self.Mx)),
                                     float(np.mean(self.My))))

    def net_M(self):
        return float(np.mean(self.Mx)), float(np.mean(self.My)), float(np.mean(self.Mz))

    def fid_arrays(self):
        if len(self._fid_sig) < 2:
            return np.array([0.0, DT]), np.array([0.0 + 0j, 0.0 + 0j])
        return np.array(self._fid_t), np.array(self._fid_sig)

    def fft_arrays(self):
        t, sig = self.fid_arrays()
        if len(sig) < 8:
            return np.array([-1.0, 0.0, 1.0]), np.array([0.0, 0.0, 0.0])
        dt  = t[1] - t[0]
        N   = len(sig)
        freq = np.fft.fftshift(fftfreq(N, d=dt))
        spec = np.fft.fftshift(np.abs(fft(sig)))
        return freq, spec

    def reset(self, N=None, dB0_ppm=None):
        if N is not None:
            self.N = N
        if dB0_ppm is not None:
            self.dB0_ppm = dB0_ppm
        self._build_offsets()
        self._init_equilibrium()


# ── 3-D Projection ────────────────────────────────────────────────────────────
def project(x, y, z, cx, cy, scale, elev, azim):
    """Simple perspective-free 3-D → 2-D projection."""
    ca, sa = math.cos(azim), math.sin(azim)
    xr =  x * ca + y * sa
    yr = -x * sa + y * ca

    ce, se = math.cos(elev), math.sin(elev)
    xs = xr
    ys = yr * se - z * ce
    zs = yr * ce + z * se
    return int(cx + xs * scale), int(cy - ys * scale), zs


def arrow3d(surf, x0, y0, z0, x1, y1, z1,
            cx, cy, scale, elev, azim,
            color, width=2, head=8):
    p0 = project(x0, y0, z0, cx, cy, scale, elev, azim)
    p1 = project(x1, y1, z1, cx, cy, scale, elev, azim)
    ax0, ay0 = p0[0], p0[1]
    ax1, ay1 = p1[0], p1[1]
    if (ax0, ay0) == (ax1, ay1):
        return
    pygame.draw.line(surf, color, (ax0, ay0), (ax1, ay1), width)
    dx, dy = ax1 - ax0, ay1 - ay0
    ln = math.hypot(dx, dy)
    if ln < 1:
        return
    ux, uy = dx / ln, dy / ln
    lx, ly = -uy, ux
    pygame.draw.polygon(surf, color, [
        (ax1, ay1),
        (int(ax1 - ux*head + lx*head*.4), int(ay1 - uy*head + ly*head*.4)),
        (int(ax1 - ux*head - lx*head*.4), int(ay1 - uy*head - ly*head*.4)),
    ])


def spin_color(mx, my):
    """Map a spin's XY direction to a hue."""
    phase = math.atan2(my, mx)
    r = int(40  + 80  * (0.5 + 0.5 * math.cos(phase)))
    g = int(160 + 60  * (0.5 + 0.5 * math.sin(phase)))
    b = int(200 - 80  * (0.5 + 0.5 * math.cos(phase)))
    return (max(0,min(255,r)), max(0,min(255,g)), max(0,min(255,b)))


# ── Drawing Routines ──────────────────────────────────────────────────────────
def draw_bloch_sphere(surf, rect, ens, elev, azim, show_individual):
    cx = rect.x + rect.w // 2
    cy = rect.y + rect.h // 2
    sc = min(rect.w, rect.h) * 0.37

    # Sphere circle
    pygame.draw.circle(surf, DIM, (cx, cy), int(sc), 1)

    # Equator ellipse
    ery = int(sc * abs(math.sin(elev)))
    if ery > 1:
        pygame.draw.ellipse(surf, DIM, (cx - int(sc), cy - ery, int(sc)*2, ery*2), 1)

    # Meridian (vertical circle in YZ plane)
    pts_m = []
    for ang in range(0, 361, 5):
        a = math.radians(ang)
        px, py, _ = project(0, math.cos(a), math.sin(a), cx, cy, sc, elev, azim)
        pts_m.append((px, py))
    if len(pts_m) > 2:
        pygame.draw.lines(surf, DIM, True, pts_m, 1)

    # Axes
    arrow3d(surf, 0,0,0, 1.3,0,0,  cx,cy,sc,elev,azim, RED_AX,   1, 7)
    arrow3d(surf, 0,0,0, 0,1.3,0,  cx,cy,sc,elev,azim, YELLOW_AX,1, 7)
    arrow3d(surf, 0,0,0, 0,0,1.4,  cx,cy,sc,elev,azim, GREEN_AX, 2, 8)

    font = pygame.font.SysFont("consolas", 11)
    for lbl, vec, col in [("x",(1.42,0,0),RED_AX),
                           ("y",(0,1.42,0),YELLOW_AX),
                           ("z",(0,0,1.55),GREEN_AX)]:
        px, py, _ = project(*vec, cx, cy, sc, elev, azim)
        surf.blit(font.render(lbl, True, col), (px-4, py-7))

    # Individual spin vectors
    if show_individual:
        for i in range(ens.N):
            mx, my, mz = ens.Mx[i], ens.My[i], ens.Mz[i]
            mag = math.sqrt(mx*mx + my*my + mz*mz)
            if mag < 0.005:
                continue
            col = spin_color(mx, my)
            arrow3d(surf, 0,0,0, mx,my,mz, cx,cy,sc,elev,azim, col, 1, 4)

    # Net M vector (amber, thick)
    nmx, nmy, nmz = ens.net_M()
    arrow3d(surf, 0,0,0, nmx,nmy,nmz, cx,cy,sc,elev,azim, ACCENT2, 3, 10)

    # Net M label
    net_mag = math.sqrt(nmx**2 + nmy**2 + nmz**2)
    px, py, _ = project(nmx, nmy, nmz, cx, cy, sc, elev, azim)
    surf.blit(font.render(f"|M|={net_mag:.2f}", True, ACCENT2), (px+4, py-8))


def draw_topdown(surf, rect, ens):
    """XY plane top-down view — dephasing is clearly visible here."""
    cx = rect.x + rect.w // 2
    cy = rect.y + rect.h // 2
    r  = int(min(rect.w, rect.h) * 0.40)

    pygame.draw.circle(surf, DIM, (cx, cy), r, 1)
    pygame.draw.line(surf, (50, 20, 20), (cx-r, cy), (cx+r, cy), 1)
    pygame.draw.line(surf, (20, 50, 20), (cx, cy-r), (cx, cy+r), 1)

    font = pygame.font.SysFont("consolas", 11)
    surf.blit(font.render("Mx →", True, RED_AX),  (cx+r-28, cy+4))
    surf.blit(font.render("↑ My", True, GREEN_AX), (cx+4,    cy-r-2))

    for i in range(ens.N):
        mx, my = ens.Mx[i], ens.My[i]
        if math.hypot(mx, my) < 0.005:
            continue
        col = spin_color(mx, my)
        ex = int(cx + mx * r)
        ey = int(cy - my * r)
        pygame.draw.line(surf, col, (cx, cy), (ex, ey), 1)
        pygame.draw.circle(surf, col, (ex, ey), 2)

    nmx, nmy, _ = ens.net_M()
    ex = int(cx + nmx * r)
    ey = int(cy - nmy * r)
    pygame.draw.line(surf, ACCENT2, (cx, cy), (ex, ey), 3)
    pygame.draw.circle(surf, ACCENT2, (ex, ey), 5)

    # Phase fan indicator
    if ens.N > 1:
        phases = [math.atan2(ens.My[i], ens.Mx[i]) for i in range(ens.N)
                  if math.hypot(ens.Mx[i], ens.My[i]) > 0.01]
        if phases:
            spread_deg = math.degrees(max(phases) - min(phases))
            surf.blit(font.render(f"spread: {spread_deg:.0f}°", True, GREY),
                      (rect.x+6, rect.y+rect.h-16))


def draw_fid(surf, rect, t_arr, fid_sig, font):
    if len(t_arr) < 2:
        return
    x0, y0, w, h = rect.x, rect.y, rect.w, rect.h
    mid = y0 + h // 2

    pygame.draw.line(surf, DIM, (x0, mid), (x0+w, mid), 1)

    sig_r = fid_sig.real
    sig_i = fid_sig.imag
    mag   = np.abs(fid_sig)
    n     = len(sig_r)
    scale = max(mag.max(), 1e-9)

    def sx_(i, v):
        px = x0 + int(i / (n-1) * w)
        py = mid - int(v / scale * h * 0.43)
        return px, py

    # Envelope shading
    env_top = [sx_(i,  mag[i]) for i in range(n)]
    env_bot = [sx_(i, -mag[i]) for i in range(n)]
    if len(env_top) >= 2:
        poly = env_top + list(reversed(env_bot))
        pygame.draw.polygon(surf, (0, 45, 38), poly)

    # Imaginary (dimmer)
    pts_i = [sx_(i, sig_i[i]) for i in range(n)]
    if len(pts_i) >= 2:
        pygame.draw.lines(surf, (0, 100, 80), False, pts_i, 1)

    # Real (bright)
    pts_r = [sx_(i, sig_r[i]) for i in range(n)]
    if len(pts_r) >= 2:
        pygame.draw.lines(surf, ACCENT, False, pts_r, 2)

    # Envelope outline
    if len(env_top) >= 2:
        pygame.draw.lines(surf, (0, 130, 100), False, env_top, 1)
        pygame.draw.lines(surf, (0, 130, 100), False, env_bot, 1)

    t_ms = t_arr[-1] * 1000
    surf.blit(font.render(f"t = {t_ms:.1f} ms", True, GREY),  (x0+5, y0+3))
    surf.blit(font.render("— Re(Mxy)",           True, ACCENT),(x0+w-80, y0+3))


def draw_fft(surf, rect, freq, spec, font):
    if len(freq) < 4:
        return
    x0, y0, w, h = rect.x, rect.y, rect.w, rect.h

    peak_i = int(np.argmax(spec))
    peak_f = freq[peak_i]

    # Zoom window around peak
    half_max = spec[peak_i] * 0.5
    above = np.where(spec >= half_max)[0]
    if len(above) >= 2:
        lw   = abs(float(freq[above[-1]]) - float(freq[above[0]]))
        zoom = max(lw * 5, 3.0)
    else:
        zoom = 15.0

    mask = (freq >= peak_f - zoom) & (freq <= peak_f + zoom)
    if not mask.any():
        mask = np.ones(len(freq), dtype=bool)

    pf   = freq[mask]
    ps   = spec[mask]
    fmin, fmax = float(pf[0]), float(pf[-1])
    smax = float(ps.max()) if ps.max() > 0 else 1.0
    fspan = max(fmax - fmin, 1e-9)

    def to_s(f, s):
        px = x0 + int((float(f) - fmin) / fspan * w)
        py = y0 + h - 4 - int(float(s) / smax * (h - 10))
        return px, py

    base = y0 + h - 4
    pts = [to_s(pf[i], ps[i]) for i in range(len(pf))]
    if len(pts) >= 2:
        poly = [(x0, base)] + pts + [(x0+w, base)]
        pygame.draw.polygon(surf, (0, 55, 45), poly)
        pygame.draw.lines(surf, ACCENT, False, pts, 2)

    # Peak marker
    pi_mask = int(np.argmax(ps))
    ppx, ppy = to_s(pf[pi_mask], ps[pi_mask])
    pygame.draw.line(surf, (255, 80, 80), (ppx, ppy-8), (ppx, ppy+8), 1)
    surf.blit(font.render(f"{peak_f:.2f} Hz", True, (255,80,80)), (ppx+3, ppy-10))

    # Axis label
    surf.blit(font.render("Freq (Hz)", True, GREY), (x0+w-68, y0+h-15))


# ── UI Widgets ────────────────────────────────────────────────────────────────
class Slider:
    def __init__(self, x, y, w, label, vmin, vmax, val, fmt="{:.2f}"):
        self.rx, self.ry = x, y + 18
        self.rw = w
        self.label = label
        self.vmin, self.vmax, self.value = vmin, vmax, val
        self.fmt = fmt
        self.dragging = False
        self._sync()

    def _sync(self):
        frac = (self.value - self.vmin) / (self.vmax - self.vmin)
        self.hx = int(self.rx + frac * self.rw)

    def draw(self, surf, font):
        surf.blit(font.render(f"{self.label}: {self.fmt.format(self.value)}",
                              True, WHITE), (self.rx, self.ry - 16))
        pygame.draw.rect(surf, SLIDER_TRK,  (self.rx, self.ry, self.rw, 6), border_radius=3)
        fw = max(0, self.hx - self.rx)
        pygame.draw.rect(surf, SLIDER_FILL, (self.rx, self.ry, fw, 6), border_radius=3)
        pygame.draw.circle(surf, WHITE, (self.hx, self.ry+3), 7)

    def handle(self, event):
        hr = pygame.Rect(self.hx-9, self.ry-5, 18, 16)
        if event.type == pygame.MOUSEBUTTONDOWN and hr.collidepoint(event.pos):
            self.dragging = True
        elif event.type == pygame.MOUSEBUTTONUP:
            self.dragging = False
        elif event.type == pygame.MOUSEMOTION and self.dragging:
            frac = max(0.0, min(1.0, (event.pos[0]-self.rx)/self.rw))
            self.value = self.vmin + frac*(self.vmax-self.vmin)
            self._sync()
            return True
        return False


class Button:
    def __init__(self, x, y, w, h, label):
        self.rect  = pygame.Rect(x, y, w, h)
        self.label = label
        self.hover = False
        self.flash = 0

    def draw(self, surf, font):
        col = BTN_ACT if self.flash > 0 else (BTN_HOV if self.hover else BTN_NRM)
        pygame.draw.rect(surf, col, self.rect, border_radius=6)
        pygame.draw.rect(surf, ACCENT, self.rect, 1, border_radius=6)
        surf.blit(font.render(self.label, True, ACCENT),
                  font.render(self.label, True, ACCENT).get_rect(center=self.rect.center))
        if self.flash > 0:
            self.flash -= 1

    def handle(self, event):
        if event.type == pygame.MOUSEMOTION:
            self.hover = self.rect.collidepoint(event.pos)
        if event.type == pygame.MOUSEBUTTONDOWN and self.rect.collidepoint(event.pos):
            self.flash = 8
            return True
        return False


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    pygame.init()
    W, H = 1300, 800
    screen = pygame.display.set_mode((W, H), pygame.RESIZABLE)
    pygame.display.set_caption("NMR Hydrogen Spin Simulation")
    clock = pygame.time.Clock()

    font_s = pygame.font.SysFont("consolas", 12)
    font_m = pygame.font.SysFont("consolas", 14)
    font_t = pygame.font.SysFont("consolas", 15, bold=True)

    CTRL_W = 270
    PAD    = 8

    def make_rects(W, H):
        vx = CTRL_W + PAD
        vw = W - vx - PAD
        vh = (H - PAD * 3) // 2
        hw = vw // 2 - PAD
        return (
            pygame.Rect(vx,       PAD,          hw, vh),
            pygame.Rect(vx+hw+PAD,PAD,          vw-hw-PAD, vh),
            pygame.Rect(vx,       PAD*2+vh,     hw, vh),
            pygame.Rect(vx+hw+PAD,PAD*2+vh,     vw-hw-PAD, vh),
        )

    bloch_r, top_r, fid_r, fft_r = make_rects(W, H)

    sx = 16
    sliders = [
        Slider(sx,  44, CTRL_W-32, "T1 (s)",       0.1,  5.0,  1.0),
        Slider(sx,  96, CTRL_W-32, "T2 (s)",       0.01, 0.5,  0.08),
        Slider(sx, 148, CTRL_W-32, "Flip angle",   1.0,  180., 90.0, "{:.0f}°"),
        Slider(sx, 200, CTRL_W-32, "Inhomog ppm",  0.0,  20.,  5.0,  "{:.1f}"),
        Slider(sx, 252, CTRL_W-32, "N spins",      4.0,  64.,  24.0, "{:.0f}"),
        Slider(sx, 304, CTRL_W-32, "Sim speed",    0.5,  10.,  3.0,  "{:.1f}x"),
    ]
    btn_pulse  = Button(sx,        360, 118, 32, "RF Pulse")
    btn_reset  = Button(sx+128,    360, 112, 32, "Reset")
    btn_toggle = Button(sx,        402, 240, 26, "Toggle individual spins")

    show_individual = True

    ens = SpinEnsemble(N=24, T1=1.0, T2=0.08, dB0_ppm=5.0, flip_deg=90.0)

    elev = math.radians(25)
    azim = math.radians(-40)
    drag_cam   = False
    drag_start = (0, 0)
    elev0 = elev
    azim0 = azim

    simulating = False

    running = True
    while running:
        clock.tick(30)

        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_ESCAPE:
                    running = False
                if event.key == pygame.K_SPACE:
                    ens.flip_deg = sliders[2].value
                    ens.T1 = sliders[0].value
                    ens.T2 = sliders[1].value
                    ens.apply_rf_pulse()
                    simulating = True

            if event.type == pygame.VIDEORESIZE:
                W, H = event.w, event.h
                screen = pygame.display.set_mode((W, H), pygame.RESIZABLE)
                bloch_r, top_r, fid_r, fft_r = make_rects(W, H)

            if event.type == pygame.MOUSEBUTTONDOWN and event.button == 3:
                if bloch_r.collidepoint(event.pos):
                    drag_cam = True
                    drag_start = event.pos
                    elev0, azim0 = elev, azim
            if event.type == pygame.MOUSEBUTTONUP and event.button == 3:
                drag_cam = False
            if event.type == pygame.MOUSEMOTION and drag_cam:
                dx = event.pos[0] - drag_start[0]
                dy = event.pos[1] - drag_start[1]
                azim = azim0 + dx * 0.012
                elev = max(-math.pi/2+0.05,
                           min(math.pi/2-0.05, elev0 - dy * 0.012))

            for sl in sliders:
                sl.handle(event)

            if btn_pulse.handle(event):
                N   = int(sliders[4].value)
                dB0 = sliders[3].value
                if N != ens.N or abs(dB0 - ens.dB0_ppm) > 0.01:
                    ens = SpinEnsemble(N=N, T1=sliders[0].value,
                                       T2=sliders[1].value,
                                       dB0_ppm=dB0,
                                       flip_deg=sliders[2].value)
                else:
                    ens.T1       = sliders[0].value
                    ens.T2       = sliders[1].value
                    ens.flip_deg = sliders[2].value
                ens.apply_rf_pulse()
                simulating = True

            if btn_reset.handle(event):
                ens = SpinEnsemble(
                    N=int(sliders[4].value),
                    T1=sliders[0].value, T2=sliders[1].value,
                    dB0_ppm=sliders[3].value, flip_deg=sliders[2].value)
                simulating = False

            if btn_toggle.handle(event):
                show_individual = not show_individual

        # ── Step simulation ────────────────────────────────────────────────
        if simulating and ens.t < T_MAX:
            speed   = sliders[5].value
            n_steps = max(1, int(speed * 4))
            ens.T1  = sliders[0].value
            ens.T2  = sliders[1].value
            for _ in range(n_steps):
                if ens.t < T_MAX:
                    ens.step(DT)
        elif simulating:
            simulating = False

        # ── Render ────────────────────────────────────────────────────────
        screen.fill(BG)
        pygame.draw.rect(screen, CTRL_BG, (0, 0, CTRL_W, H))

        screen.blit(font_t.render("H Spin Simulation", True, ACCENT), (sx, 14))

        for sl in sliders:
            sl.draw(screen, font_s)
        btn_pulse.draw(screen, font_m)
        btn_reset.draw(screen, font_m)
        btn_toggle.draw(screen, font_s)

        # Status
        nmx, nmy, nmz = ens.net_M()
        nxy = math.hypot(nmx, nmy)
        status = [
            f"t = {ens.t*1000:.1f} ms",
            f"|Mxy| = {nxy:.3f}",
            f"Mz    = {nmz:.3f}",
            f"{'▶ simulating' if simulating else '■ done / paused'}",
            "",
            "SPACE or RF Pulse = tip spins",
            "Right-drag sphere = rotate",
            "Reset = back to equilibrium",
            "Toggle = show/hide spins",
        ]
        for i, line in enumerate(status):
            col = ACCENT if (i == 3 and simulating) else (GREY if i > 4 else WHITE)
            screen.blit(font_s.render(line, True, col), (sx, 442 + i*17))

        # Panels
        for r in [bloch_r, top_r, fid_r, fft_r]:
            pygame.draw.rect(screen, PANEL_BG, r, border_radius=4)

        draw_bloch_sphere(screen, bloch_r, ens, elev, azim, show_individual)
        draw_topdown(screen, top_r, ens)

        t_arr, fid_sig = ens.fid_arrays()
        draw_fid(screen, fid_r, t_arr, fid_sig, font_s)

        freq, spec = ens.fft_arrays()
        draw_fft(screen, fft_r, freq, spec, font_s)

        panel_info = [
            (bloch_r, "Bloch Sphere  (right-drag = rotate)"),
            (top_r,   "XY Precession  (top view)"),
            (fid_r,   "FID Signal"),
            (fft_r,   "FFT Spectrum"),
        ]
        for r, lbl in panel_info:
            pygame.draw.rect(screen, (10, 10, 28), (r.x, r.y, r.w, 20))
            screen.blit(font_s.render(lbl, True, ACCENT), (r.x+6, r.y+3))
            pygame.draw.rect(screen, ACCENT, r, 1, border_radius=4)

        pygame.display.flip()

    pygame.quit()


if __name__ == "__main__":
    main()