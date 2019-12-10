GlowScript
2.9
VPython

win = 500
CC = 0
G = 6.67e-11

planet_side = 3
Natoms = 2 * planet_side ** 3

L = 10e7
gray = color.gray(0.7)
mass = 1e23
Ratom = 1.5e6
k = 1.4E-23
T = 300
dt = 5
density = mass / (4 / 3 * pi * Ratom ** 3)
print(density)

center1 = vec(-L / 4, -L / 4, -L / 4)
center2 = vec(L / 4, L / 4, L / 4)

animation = canvas(width=win, height=win, align='left')
animation.range = L
animation.title = 'A "hard-sphere" gas'
s = """  Theoretical and averaged speed distributions (meters/sec).
  Initially all atoms have the same speed, but collisions
  change the speeds of the colliding atoms. One of the atoms is
  marked and leaves a trail so you can follow its path.

"""
energy_d = gdisplay(xtitle='time', ytitle='kinetic energy', align='right')
energy_plot = gcurve(gdisplay=energy_d)
animation.caption = s

d = L / 2 + Ratom
r = 0.005

Atoms = []
p = []
apos = []
pavg = mass * 1e-3  # average kinetic energy p**2/(2mass) = (3/2)kT

gmax = mass ** 2 * G / (2 * r) ** 2

for j1 in range(planet_side):
    for j2 in range(planet_side):
        for j3 in range(planet_side):
            posit = center1 + 3 * Ratom * vec(j1, j2, j3)
            Atoms.append(sphere(pos=posit, radius=Ratom, color=gray))
            apos.append(posit)
for j1 in range(planet_side):
    for j2 in range(planet_side):
        for j3 in range(planet_side):
            posit = center2 + 3 * Ratom * vec(j1, j2, j3)
            Atoms.append(sphere(pos=posit, radius=Ratom, color=gray))
            apos.append(posit)

for i in range(Natoms):
    p.append(vector(0, 0, 0))

deltav = 100  # binning for v histogram


def barx(v):
    return int(v / deltav)  # index into bars array


nhisto = int(4500 / deltav)
histo = []
for i in range(nhisto): histo.append(0.0)
histo[barx(pavg / mass)] = Natoms

gg = graph(width=win, height=0.4 * win, xmax=3000, align='left',
           xtitle='speed, m/s', ytitle='Number of atoms', ymax=Natoms * deltav / 1000)

theory = gcurve(color=color.cyan)
dv = 10
for v in range(0, 3001 + dv, dv):  # theoretical prediction
    theory.plot(v, (deltav / dv) * Natoms * 4 * pi * ((mass / (2 * pi * k * T)) ** 1.5) * exp(
        -0.5 * mass * (v ** 2) / (k * T)) * (v ** 2) * dv)

accum = []
for i in range(int(3000 / deltav)): accum.append([deltav * (i + .5), 0])
vdist = gvbars(color=color.red, delta=deltav)


def interchange(v1, v2):  # remove from v1 bar, add to v2 bar
    barx1 = barx(v1)
    barx2 = barx(v2)
    if barx1 == barx2:  return
    if barx1 >= len(histo) or barx2 >= len(histo): return
    histo[barx1] -= 1
    histo[barx2] += 1


def checkCollisions():
    hitlist = []
    r2 = 2 * Ratom
    r2 *= r2
    for i in range(Natoms):
        ai = apos[i]
        for j in range(i):
            aj = apos[j]
            dr = ai - aj
            if mag2(dr) < r2: hitlist.append([i, j])
    return hitlist


nhisto = 0  # number of histogram snapshots to average


def energyplot():
    e = 0
    for mom in p:
        e += .5 * mass * mag2(mom / mass)
    return e


gravmag_disp = gdisplay(xtitle='time', ytitle='magnitude of gravity forces total', align='right')
gravmag_plot = gcurve(gdisplay=gravmag_disp)
counter = 0

while True:
    rate(100000000)
    counter += 1
    energy_plot.plot(dt * counter, energyplot())

    # Accumulate and average histogram snapshots
    for i in range(len(accum)): accum[i][1] = (nhisto * accum[i][1] + histo[i]) / (nhisto + 1)
    if nhisto % 10 == 0:
        vdist.data = accum
    nhisto += 1

    # Update all positions
    for i in range(Natoms): Atoms[i].pos = apos[i] = apos[i] + (p[i] / mass) * dt
    gravmag = 0
    for i in range(Natoms):
        for j in range(Natoms):
            if i == j:
                continue
            atomi = apos[i]
            atomj = apos[j]
            if mag2(atomi - atomj) == 0:
                print('oops')
                print(i)
                print(j)
                break
            f = G * mass ** 2 / mag2(atomi - atomj) * norm(atomi - atomj)
            if mag(f) > gmax:
                f = norm(f) * gmax
            gravmag += mag(f)
            p[j] += f * dt
    gravmag_plot.plot(dt * counter, gravmag)

    # Check for collisions
    hitlist = checkCollisions()

    # If any collisions took place, update momenta of the two atoms
    for ij in hitlist:
        i = ij[0]
        j = ij[1]
        ptot = p[i] + p[j]
        posi = apos[i]
        posj = apos[j]
        vi = p[i] / mass
        vj = p[j] / mass
        vrel = vj - vi

        a = vrel.mag2
        if a == 0: continue;  # exactly same velocities
        rrel = posi - posj
        if rrel.mag > Ratom:
            continue  # one atom went all the way through another

        # theta is the angle between vrel and rrel:
        dx = dot(rrel, vrel.hat)  # rrel.mag*cos(theta)
        dy = cross(rrel, vrel.hat).mag  # rrel.mag*sin(theta)
        # alpha is the angle of the triangle composed of rrel, path of atom j, and a line
        #   from the center of atom i to the center of atom j where atome j hits atom i:
        alpha = asin(dy / (2 * Ratom))
        d = (2 * Ratom) * cos(alpha) - dx  # distance traveled into the atom from first contact
        deltat = d / vrel.mag  # time spent moving from first contact to position inside atom
        posi = posi - vi * deltat  # back up to contact configuration
        posj = posj - vj * deltat
        mtot = 2 * mass
        pcmi = p[i] - ptot * mass / mtot  # transform momenta to cm frame
        pcmj = p[j] - ptot * mass / mtot
        rrel = norm(rrel)
        pcmi = (pcmi - 2 * pcmi.dot(rrel) * rrel) * CC  # bounce in cm frame
        pcmj = (pcmj - 2 * pcmj.dot(rrel) * rrel) * CC
        p[i] = pcmi + ptot * mass / mtot  # transform momenta back to lab frame
        p[j] = pcmj + ptot * mass / mtot
        apos[i] = posi + (p[i] / mass) * deltat  # move forward deltat in time
        apos[j] = posj + (p[j] / mass) * deltat
        interchange(vi.mag, p[i].mag / mass)
        interchange(vj.mag, p[j].mag / mass)



