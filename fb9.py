import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pylab as mp
from matplotlib.animation import FuncAnimation
from scipy.fft import fft, ifft, fftshift
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from scipy import optimize
import math
import time
import xlsxwriter

def fftfreqlmd(N,T):
    t = np.linspace(-T / 2, T / 2, N)
    dt = t[1] - t[0]
    w1 = 2 * np.pi * np.arange(0, N, 1) / T
    w12 = fftshift(w1)
    w = w12 - w12[0]
    #w3 = fftfreq(n1, twidth) * 2 * np.pi * N
    lmd2 = 2 * np.pi * c / (-w + w0) / 1e-9  # wavelength
    lmd = fftshift(lmd2)
    return [t, dt, w, lmd]
def passive(y0, pfparam, ig):

    Ld = pfparam[-1]
    bet2, gam = pfparam[1:-1]
    elname = pfparam[0]
    Disper = 1j * (bet2[0] * (w ** 2) / 2 + bet2[1] * (w ** 3) / 6)
    # + bet2[2] * (w ** 4) / 24 + bet2[3] * (w ** 5) / 120 + bet2[4] * (w ** 6) / 720 + bet2[5] * (w ** 7) / 5040)

    def pf(y0, z):  # passive fiber

        y = y0[0:n1] + 1j * y0[n1:]
        Ax = ifft(y)
        Gw = fft((np.abs(Ax) ** 2) * Ax)
        dy1 = Disper * y + gam * 1j * Gw
        dy = np.hstack([np.real(dy1), np.imag(dy1)])
        return dy

    y02 = np.real(y0)
    y02 = np.hstack([y02, np.imag(y0)])
    sol = odeint(pf, y02, Ld, args=(), atol=1e-6, rtol=1e-8, hmax = 0.02)
    y2 = sol[-1,:]
    y1 = y2[0:n1] + 1j * y2[n1:]
    solgz = sol[:, 0:n1] + 1j * sol[:, n1:]
    ygz = ifft(solgz)
    #Ez = np.sum(np.abs(ygz)**2, 1) * dt * 1e-12 / Psg / Tr
    Ez0 = np.sum(np.abs(ygz) ** 2, 1) * dt * 1e-12
    Pz = np.abs(ygz) ** 2
    Pz0 = np.max(Pz, 1)
    return [y1, [elname, Ld, Pz0], [elname, Ld, Ez0]]
def afgain(OOm2):
    gain1 = (1 - (w / OOm2) **2)
    gain1[gain1 < 0] = 0
    return gain1
def active(y0, afparam, ig):

    #y0 - pulse spectrum at active fiber entrance
    La = afparam[-1] # active fiber length
    # active:[beta, gama1, OOm2, g0, Ez, Psg, Tr, direction, La1],
    bet2, gam, OOm2, g0, Ez, Psg, Tr, drctn = afparam[1:-1]
    elname = afparam[0]
    if elname == 'activefB':
        #g0 = getnalmg(g00, ig)
        print('nalmg0=', g0)
    #else:
        #print('resonator')
        #g0 = getgf(g00, ig)

    # bet2 - dispertion, gam - nonlinearity, OOm2 - amplif bandwidth coeff, g0 - small signal gain, EzXt, Ez - pulse energy divided on Psg and Tr, Psg - active fiber power of saturation, Tr - round trip time, drctn - direction of moving (1 - forward, 2 - back)
    Disper = 1j * (bet2[0] * (w ** 2) / 2 + bet2[1] * (w ** 3) / 6 + bet2[2] * (w ** 4) / 24 + bet2[3] * (w ** 5) / 120 + bet2[4] * (w ** 6) / 720 + bet2[5] * (w ** 7) / 5040)
    gain = afgain(OOm2)
    #print(g0)

    if drctn == -1:
        Ez = Ez[-1::-1] # reverse pulse energy data through active fiber

    x2 = np.linspace(La[0], La[-1], 4)
    uf = interp1d(La, Ez, kind='cubic')
    uf2 = uf(x2)
    Xt = tcf(x2, uf2) # Taylor coeff for Ez profile through active fiber length

    def af(y0, z):  # active fiber
        #print([z, La[0], La[-1]])
        y = y0[0:n1] + 1j * y0[n1:] # solver can't accept complex data. So y0 is real, but y is complex. The y is pulse spectrum
        Ax = ifft(y) # Pulse profile
        Gw = fft((np.abs(Ax) ** 2) * Ax) # fft from Shredinger equation nonlinear part
        Egz0 = np.sum(np.abs(Ax) ** 2) * dt * 1e-12 / Psg / Tr
        Egz1 = tvl(Xt, z, Ez[0], La[0]) / Psg / Tr
        Egz = Egz0 + Egz1
        gt = g0 / (1 + Egz) # saturated gain
        gUwx = fft(gt * Ax)
        dy1 = Disper * y + gam * 1j * Gw + gain * gUwx
        dy = np.hstack([np.real(dy1), np.imag(dy1)])
        return dy

    y02 = np.real(y0)
    y02 = np.hstack([y02, np.imag(y0)])

    sol = odeint(af, y02, La, args=(), atol=1e-10, rtol=1e-8, hmax=0.005)
    y2 = sol[-1, :]
    y1 = y2[0:n1] + 1j * y2[n1:]
    solgz = sol[:, 0:n1] + 1j * sol[:, n1:]
    ygz = ifft(solgz)
    #print(3)
    #print(np.abs(ygz))
    #Ez = np.sum(np.abs(ygz)**2, 1) * dt * 1e-12 / Psg / Tr
    Ez0 = np.sum(np.abs(ygz) ** 2, 1) * dt * 1e-12
    Pz = np.abs(ygz) ** 2
    Pz0 = np.max(Pz, 1)

    if drctn == -1:
        Ez0 = Ez0[-1::-1]
        Pz0 = Pz0[-1::-1]

    return [y1, [elname, La, Pz0], [elname, La, Ez0]]
def sesam(yw, sparam, ig):
    q0, Pa, ta, = sparam[1:]
    elname = sparam[0]
    ui = ifft(yw)
    uf = interp1d(t, np.abs(ui), kind='cubic')

    def sesamt(q, ts):
        uii = uf(ts)
        dq = -1 * (q - q0) / ta - 1 * q * np.abs(uii) ** 2 / Pa / ta
        return dq


    if ta == 0:
        q = q0 / (1 + np.abs(ui[0:-1]) ** 2 / Pa)
    else:
        q = odeint(sesamt, sparam[0], t[0:-1], atol=1e-6, rtol=1e-6, hmax=0.001)
        q = q.T
        q = q[0]

    q = np.hstack([q, q[-1]])
    yp = ((1 - q) ** 0.5) * ui
    y = fft(yp)

    E1 = np.sum(np.abs(ui) ** 2) * dt * 1e-12
    E2 = np.sum(np.abs(yp) ** 2) * dt * 1e-12
    P1 = np.max(np.abs(ui) ** 2)
    P2 = np.max(np.abs(yp) ** 2)

    return [y, [elname, np.array([0, 5]), np.array([P1, P2])], [elname, np.array([0, 5]), np.array([E1, E2])], q]
def coupler(yw, cparam, ig):
    ratio = cparam[-1]
    elname = cparam[0]
    y = ratio**0.5 * yw
    ui = ifft(yw)
    yp = ifft(y)

    E1 = np.sum(np.abs(ui) ** 2) * dt * 1e-12
    E2 = np.sum(np.abs(yp) ** 2) * dt * 1e-12
    P1 = np.max(np.abs(ui) ** 2)
    P2 = np.max(np.abs(yp) ** 2)

    return [y, [elname, np.array([0, 5]), np.array([P1, P2])], [elname, np.array([0, 5]), np.array([E1, E2])]]
def filterd(par):
    lcenter = par[0]
    fwidth = par[1]
    wcenter = w0 - 2 * np.pi * c / lcenter / 1e-9
    fltrw = np.exp(-((w - wcenter) ** 2 / 2 / fwidth ** 2) ** 10)
    return fltrw
def filter(y0, fparam, ig):
    fpr = fparam[1:]
    elname = fparam[0]
    y = filterd(fpr) * y0

    ui = ifft(y0)
    yp = ifft(y)

    E1 = np.sum(np.abs(ui) ** 2) * dt * 1e-12
    E2 = np.sum(np.abs(yp) ** 2) * dt * 1e-12
    P1 = np.max(np.abs(ui) ** 2)
    P2 = np.max(np.abs(yp) ** 2)

    return [y, [elname, np.array([0, 5]), np.array([P1, P2])], [elname, np.array([0, 5]), np.array([E1, E2])]]
def nalm(yw, nalmparam, ig):
    global yfnalm, ybnalm

    elname = nalmparam[0]
    yw01 = (0.6**0.5) *yw
    yw02 = 1j*(0.4**0.5) *yw
    passiveNALM1, activeNALM1, passiveNALM2 = nalmparam[1:]
    print('Lnalm2 = ', passiveNALM2[-1][-1])

    # nalm direction +
    y11 = passiveNALM1[0](yw01, passiveNALM1[1:], ig)
    y = y11[0]
    activeNALM1[9] = 1
    y11 = activeNALM1[0](y, activeNALM1[1:], ig)
    y = y11[0]
    activeNALM1[6] = y11[2][2]
    Ef1 = y11[2][2]
    y11 = passiveNALM2[0](y, passiveNALM2[1:], ig)
    yf = y11[0]

    # nalm direction -
    y11 = passiveNALM2[0](yw02, passiveNALM2[1:], ig)
    y = y11[0]
    activeNALM1[9] = -1
    y11 = activeNALM1[0](y, activeNALM1[1:], ig)
    y = y11[0]
    activeNALM1[6] = y11[2][2]
    Efb = y11[2][2]
    y11 = passiveNALM1[0](y, passiveNALM1[1:], ig)
    yb = y11[0]

    yfnalm = ifft(yf)*(0.6**0.5)
    ybnalm = ifft(yb)*1j*(0.4**0.5)
    Enalm = (Ef1 + Efb)/np.max(Ef1 + Efb)

    #activeNALM1[6] = (Enalm1 + Enalm2[-1::-1])/2
    ui = ifft(yw)
    yp = (ifft(yf)*(0.6**0.5) + 1j * ifft(yb)*(0.4**0.5))

    E1 = np.sum(np.abs(ui) ** 2) * dt * 1e-12
    E2 = np.sum(np.abs(yp) ** 2) * dt * 1e-12
    P1 = np.max(np.abs(ui) ** 2)
    P2 = np.max(np.abs(yp) ** 2)
    y = fft(yp)

    return [y, [elname, np.array([0, 5]), np.array([P1, P2])], [elname, np.array([0, 5]), np.array([E1, E2])], [yfnalm, ybnalm]]
def cavitystruct():

    Aaf = np.pi * (9 *1e-6) ** 2 / 4 # square of fiber core
    n2 = 2.5e-20
    gama1 = n2 * w0 / c / Aaf

    # active fiber
    beta = np.array([0.069551743234435, 0.000288051116037, -0.000001103550617, 0.000000005649610, -0.000000000034708, 0.000000000000249])
    OOm2 = 20 # gain width
    g01 = 20
    g02 = 20
    Psg = 10e-3
    Tr = 1.3/2e8
    La1 = np.linspace(0, 1, 500)
    La2 = np.linspace(0, 1, 500)
    Ez = 1 * 1000 * np.ones(len(La1))
    direction = 1
    # end

    # passive fiber
    betp = np.array([6.276e-2, 8.119e-5, -1.321e-7, -3.032e-10, -4.196e-13, -2.570e-16])
    betp = beta
    Aaf = np.pi * (9 *1e-6) ** 2 / 4 # square of fiber core
    n2 = 2.5e-20
    gamp1 = n2 * w0 / c / Aaf
    Lp1 = np.linspace(0, 150, 1000)
    Lp2 = np.linspace(0, 170, 1000)
    Lp3 = np.linspace(0, 2.2, 100)
    Lp4 = np.linspace(0, 25, 1000)
    # coupler
    cplrratio = 0.5

    # sesam
    q0 = 0.5
    Pa = 20
    ta = 0

    # filter
    lcenter = lmc # filter center in wavelength
    fwidth = 3 # filter width
    # end filter
    """
    cavelm = {passive: [betp, gamp1, Lp1],
              active: [beta, gama1, OOm2, g0, Ez, Psg, Tr, direction, La1],
              coupler: cplrratio,
              sesam: [q0, Pa, ta],
              filter: [lcenter, fwidth]}

    cavstr = [filter, passive, coupler, active, filter, sesam]
    """
    passiveA = [passive, 'passivefA', betp, gamp1, Lp1]
    passiveB = [passive, 'passivefB', betp, gamp1, Lp2]
    passiveC = [passive, 'passivefC', betp, gamp1, Lp3]
    passiveD = [passive, 'passivefD', betp, gamp1, Lp4]
    activeA = [active, 'activefA', beta, gama1, OOm2, g01, Ez, Psg, Tr, direction, La1]
    activeB = [active, 'activefB', beta, gama1, OOm2, g02, Ez, Psg, Tr, direction, La2]
    filterB = [filter, 'filterB', lcenter, fwidth]
    couplerA = [coupler, 'couplerA', cplrratio]
    sesamA = [sesam, 'sesamA', q0, Pa, ta]
    nalmA = [nalm, 'nalmA', passiveC, activeB, passiveD]

    cavstruct = [activeA, passiveA, nalmA, couplerA, passiveB, filterB]
    someparams = [Lp1[-1], OOm2, [lcenter, fwidth], g01]

    return cavstruct, someparams
def initdata(twidth, n1):



    c = 3 * 1e8 / 1e12 # speed of light in m/ps
    lmc = 1930 # active fiber amplification band center
    t, dt, w, lmd, = fftfreqlmd(n1, twidth,) # define t, w, lmd scale. lmd - wavelength scale

    y = 0.3*np.exp(-(t / 100) ** 2) # seed pulse
    yf = fft(y) # seed spectrum
    yff = fftshift(yf)

    return (yf, t, dt, w, lmd)
def initxlsfile():
    # (lmd, yspectr, t, [ytime, ytcomp], tchirp, ychirp, [TBPt, TBPtcomp], TBPlmd, [TBPdata, TBPdatacomp], [TBPt60, TBPtcomp60], TBPlmd60, Ip, Ep)
    worksheet.set_column('B:L', 20)
    bold = workbook.add_format({'bold': True})
    worksheet.write('A1', 'Number', bold)
    worksheet.write('B1', '-3dB durtn, ps', bold)
    worksheet.write('C1', 'compr -3dB durtn, ps', bold)
    worksheet.write('D1', '-3dB spctrm wdth, nm', bold)
    worksheet.write('E1', 'TBP', bold)
    worksheet.write('F1', 'compr TBP', bold)
    worksheet.write('G1', '-60dB durtn, ps', bold)
    worksheet.write('H1', 'compr -60dB durtn, ps', bold)
    worksheet.write('I1', '-60dB spctrm wdth, nm', bold)
    worksheet.write('J1', 'peak intst, W', bold)
    worksheet.write('K1', 'enrj, nJ', bold)
    worksheet.write('L1', 'nalm passiveD L, m', bold)

    return 0
def initgraph():
    axs[0].set_xlim([lmd[0], lmd[-1]])  # axs[0] - to show spectrum
    axs[0].set_ylim([-100, 0.1])
    axs[1].set_xlim([t[0], t[-1]])  # axs[1] - to show pulse profile in time
    axs[1].set_ylim([-100, 0.1])

    axs[2].set_ylim(0, 1)
    axs[2].set_title('Peak intensity, W')
    axs[2].set_xlabel('Cavity round trip')
    axs[2].set_xlim(0, frames0)  # axs[2] - to show pulse parameter (peak intensity for exmp.) with cavity round trip

    axs[3].set_xlim(-20,
                    400)  # axs[3] - to show some cavity element or pulse param as pulse goes through cavity element. now its g (saturated ampligication) through active fiber length
    axs[3].set_ylim(0, 50)
    axs[3].set_title('Peak intensity, W')
    axs[3].set_xlabel('Cavity coordinate, m')

    axch.set_ylim([1900, 2000])

    axsnalm.set_xlim([t[0], t[-1]])
    axsnalm.set_ylim([-3, 3])
    axsnalm.set_title('Pulse phaze difference at NALM output, pi multiple')
    axsnalm.set_xlabel('Pulse time delay, ps')

    #axsnalmE.set_xlim([0, 1])
    #axsnalmE.set_ylim([0, 1])
    #axsnalmE.set_title('E pulse at NALM output, rad')
    #axsnalmE.set_xlabel('NALM active fber length, m')

    line1nalm, = axsnalm.plot([], []) # 33 for nalm faze calc
    #line2nalm, = axsnalm.plot([], []) # 34 for nalm faze calc
    #line1nalmEp, = axsnalmE.plot([], [])  # 34 for nalm Ep

    line11, = axs[0].plot([], [])  # 0 spectr
    line12, = axs[0].plot([], [])  # 1 gain
    line13, = axs[0].plot([], [])  # 2 filter
    line14, = axs[0].plot([], [], 'o')  # 3 TBPpointsw
    line21, = axs[1].plot([], [])  # 4 profile in time domain
    line22, = axch.plot([], [], 'r')  # 5 chirp
    line23, = axs[1].plot([], [], 'o')  # 6 TBPpointsp
    line24, = axs[1].plot([], [])  # 7 profile compressed
    line25, = axs[1].plot([], [], 'o')  # 8 profile compressed TBP points
    line26, = axs[1].plot([], [])  # 9 ...
    line27, = axs[1].plot([], [])  # 10 ...
    line31, = axs[2].plot([], [])  # 11 pulse peak power at cavity round times
    line32, = axs[2].plot([], [])  # 12 ...
    line33, = axs[2].plot([], [])  # 13 ...

    line41, = axs[3].plot([], [])  # 14 passive fiber energy
    line42, = axs[3].plot([], [])  # 15 coupler energy
    line43, = axs[3].plot([], [])  # 16 active fiber energy
    line44, = axs[3].plot([], [])  # 17 filter energy
    line45, = axs[3].plot([], [])  # 18 sesam energy
    line46, = axs[3].plot([], [])  # 19
    line47, = axs[3].plot([], [])  # 20
    line48, = axs[3].plot([], [])  # 21
    line49, = axs[3].plot([], [])  # 22
    line411, = axs[3].plot([], [])  # 23
    line412, = axs[3].plot([], [])  # 24

    line413, = axs[3].plot([], [])  # 25 passive fiber peak int
    line414, = axs[3].plot([], [])  # 26 coupler peak int
    line415, = axs[3].plot([], [])  # 27 active fiber peak int
    line416, = axs[3].plot([], [])  # 28 filter peak int
    line417, = axs[3].plot([], [])  # 29 sesam peak int
    line418, = axs[3].plot([], [])  # 30
    line419, = axs[3].plot([], [])  # 31
    line4111, = axs[3].plot([], [])  # 32


    axs[0].legend(['Spectrum', 'Gain', 'Filter', '-3dB level'], loc='lower left')
    axs[1].legend(['Intensity', '-3dB Intsty', 'Int compressed', '-3dB Int cmprsd'], loc='lower left')
    axch.legend(['Chirp, nm'], loc='upper right')
    leg3 = []
    for ity in cvtstruct:
        leg3.append(ity[1])

    #axs[3].legend(['psv fiber', 'coupler', 'act fiber', 'filter', 'sesam'])
    axs[3].legend(leg3)

    return (line11, line12, line13, line14, line21, line22, line23, line24, line25,
            line26, line27, line31, line32, line33, line41, line42, line43, line44,
            line45, line46, line47, line48, line49, line411, line412, line413, line414,
            line415, line416, line417, line418, line419, line4111,
            line1nalm,)
def FWHM3dB2(t, y): # get pulse width at -3 dB
    y = np.abs(y)
    uf = interp1d(t, y, kind='cubic')
    iy = []
    iy1 = []
    y2 = np.abs(y) - 0.5 * np.max(np.abs(y))
    for ix in range(0, len(np.abs(y2))):
        if y2[ix] > 0:
            iy.append(ix)
    t3 = [t[iy[0] - 1], t[iy[0]]]
    yp3 = [y[iy[0] - 1], y[iy[0]]]
    k1 = (yp3[-1] - yp3[0]) / (t3[-1] - t3[0])
    b1 = yp3[0] - t3[0] * k1
    t3x = (0.5 * np.max(np.abs(y)) - b1) / k1
    y3x = 0.5 * np.max(np.abs(y))

    t4 = [t[iy[-1]], t[iy[-1]+1]]
    yp4 = [y[iy[-1]], y[iy[-1]+1]]
    k1 = (yp4[-1] - yp4[0]) / (t4[-1] - t4[0])
    b1 = yp4[0] - t4[0] * k1
    t4x = (0.5 * np.max(np.abs(y)) - b1) / k1
    y4x = 0.5 * np.max(np.abs(y))
    tx = [t3x, t4x]
    yx = [uf(t3x), uf(t4x)]
    return np.array([tx, yx])  # for del
def FWHM3dB(t, y): # get pulse width at -3 dB
    y = np.abs(y)
    ti = np.linspace(t[0], t[-1], len(t)*100)
    uf = interp1d(t, y, kind='cubic')
    iy = []
    ufi = uf(ti)
    y2 = np.abs(ufi) - 0.5 * np.max(np.abs(ufi))
    for ix in range(0, len(np.abs(y2))):
        if y2[ix] > 0:
            iy.append(ix)
    txi = ti[iy]
    if len(txi) >= 2:
        tx = [ti[iy[0]], ti[iy[-1]]]
        yx = [ufi[iy[0]], ufi[iy[-1]]]
    else:
        tx = np.array([0, 1])
        yx = np.array([1, 1])
    return np.array([tx, yx])
def FWHM60dB(t, y): # get pulse width at -3 dB
    y = np.abs(y)
    ti = np.linspace(t[0], t[-1], len(t)*100)
    uf = interp1d(t, y, kind='cubic')
    iy = []
    ufi = uf(ti)
    y2 = 10 * np.log10(np.abs(ufi)/np.max(np.abs(ufi))) + 60
    for ix in range(0, len(np.abs(y2))):
        if y2[ix] > 0:
            iy.append(ix)
    txi = ti[iy]
    if len(txi) >= 2:
        tx = [ti[iy[0]], ti[iy[-1]]]
        yx = [ufi[iy[0]], ufi[iy[-1]]]
    else:
        tx = np.array([0, 1])
        yx = np.array([1, 1])
    return np.array([tx, yx])
def tbp(t, w, yw):
    y = ifft(yw)
    yv = fftshift(yw)
    tp = FWHM3dB(t, y ** 2)[0]
    wp = FWHM3dB(fftshift(w), yv ** 2)[0]
    dt = tp[-1] - tp[0]
    dw = wp[-1] - wp[0]
    tbp = dt * dw / 2 / np.pi
    return tbp
def pulseparams(ywtlmd):

    #  ywtlmd = (y, t, lmd, w)
    #w = ywtlmd[3]
    #lmd = ywtlmd[2]
    #t = ywtlmd[1]
    yt = ifft(ywtlmd)
    yw0 = fftshift(ywtlmd)
    yspectr = 10 * np.log10((np.abs(yw0)**2) / np.max(np.abs(yw0)**2)) # spectrum intensity
    ytime = 10 * np.log10((np.abs(yt)**2) / np.max(np.abs(yt)**2)) # profile intensity

    itchirp = np.array([ix for ix in range(len(ytime)) if ytime[ix] > -60])
    #print("len(chirp) at 268 line=", len(itchirp))
    if len(itchirp)>1:
        tchirp = t[itchirp]
        ii = min(np.abs(t))
        tlist = t.tolist()
        it = tlist.index(ii)
        ch1 = np.unwrap(np.angle(yt[it:]))
        ch2 = np.unwrap(np.angle(np.flip(yt[:it])))
        ch = np.hstack([np.flip(ch2),ch1])

        chirp0 = -1 * np.diff(ch) / (t[1]-t[0])
        chirp0 = np.hstack([chirp0, chirp0[-1]])
        chirp = chirp0[itchirp]
        ychirp = 2 * np.pi * c / (w0 + chirp) /(1e-9)
    else:
        tchirp = 0
        ychirp = 0

    #def tbp(t, w, yw):
    TBPdata = tbp(t, w, ywtlmd) #def tbp(t, w, yw):
    TBPt = FWHM3dB(t, ifft(ywtlmd)**2)[0] #
    TBPt60 = FWHM60dB(t, ifft(ywtlmd) ** 2)[0]  #
    TBPw = FWHM3dB(fftshift(w), np.abs(fftshift(ywtlmd))**2)[0]
    TBPw60 = FWHM60dB(fftshift(w), np.abs(fftshift(ywtlmd)) ** 2)[0]
    TBPlmd = 2 * np.pi * c / (-TBPw + w0) / 1e-9
    TBPlmd60 = 2 * np.pi * c / (-TBPw60 + w0) / 1e-9
    disptocompress = optimize.fmin(compresspulse,  disptocompress0[0], args = (ywtlmd,) ,disp = False)
    disptocompress0[0] = disptocompress

    ywcomp0 = getpulsecompressed(ywtlmd, disptocompress)
    ywcomp0 = movetotime(ywcomp0, 0)
    ytcomp0 = np.abs(ifft(ywcomp0))
    ytcomp = 10 * np.log10((ytcomp0 ** 2) / np.max(ytcomp0 ** 2))  # profile intensity

    TBPdatacomp = tbp(t, w, ywcomp0) #def tbp(t, w, yw):
    TBPtcomp = FWHM3dB(t, ifft(ywcomp0)**2)[0] #
    TBPtcomp60 = FWHM60dB(t, ifft(ywcomp0) ** 2)[0]  #

    Ep = np.sum(np.abs(yt) ** 2) * dt * 1e-12
    Ip = np.max(np.abs(yt) ** 2)

    print("time-bandwidth product", TBPdata)
    print("time-bandwidth product comp", TBPdatacomp)
    print('-3dB pulse duration, ps=', np.abs(TBPt[1]-TBPt[0]))
    print('-3dB pulse duration comp, ps=', np.abs(TBPtcomp[1] - TBPtcomp[0]))
    print('-3dB spectrum width, nm=', np.abs(TBPlmd[1] - TBPlmd[0]))

    pp = (lmd, yspectr, t, [ytime, ytcomp], tchirp, ychirp, [TBPt, TBPtcomp], TBPlmd, [TBPdata, TBPdatacomp], [TBPt60, TBPtcomp60], TBPlmd60, Ip, Ep)

    # tp, yp, lmdp, ywp, tbp - time points, lmd points for tbp
#    print(tchirp)
#    print(ychirp)
    return pp
def report(ywtlmd, Pp, Ip, Ep):
    # pdata = (lmd, yspectr, t, [ytime, ytcomp], tchirp, ychirp, [TBPt, TBPtcomp], TBPlmd, [TBPdata, TBPdatacomp])
    # cavityparams = (gain,)

    global yfnalm, ybnalm

    pdata = pulseparams(ywtlmd)  # to get pdata
    gain = afgain(somepar[1])
    gain = 10 * np.log10((np.abs(gain) ** 2) / np.max(np.abs(gain) ** 2))
    gain = fftshift(gain)
    # end gain

    # get filter
    flt = filterd(somepar[2])
    flt = 10 * np.log10((np.abs(flt) ** 2) / np.max(np.abs(flt) ** 2))
    flt = fftshift(flt)
    # end filter

    axch.set_ylim([np.min(pdata[5]), np.max(pdata[5])])
    #axs[3].set_xlim([0, np.max(np.array(Ip[0][1])*1.5)])
    #axs[3].set_ylim([0, np.max(np.array(Ip[2][2])*1.5)])
    #axs[2].set_xlim([0, npp[-1]+1])
    #axs[2].set_ylim([0, np.max(Pp)])
    npp = np.arange(1, len(Pp) + 1)

    TBPlmd = pdata[7]
    TBPt = pdata[6][0]
    TBPtcomp = pdata[6][1]

    lines[0].set_data(pdata[0], pdata[1])  # 0 spectrum intesity log10
    lines[1].set_data(pdata[0], gain)  # gain 1
    lines[2].set_data(pdata[0], flt)  # filter 1
    lines[3].set_data(TBPlmd, np.array([-3,-3]))  # TBPpoints lmd

    lines[4].set_data(pdata[2], pdata[3][0])  # 4 profile intencity log10
    lines[5].set_data(pdata[4], pdata[5])  # 5 chirp
    lines[6].set_data(TBPt, np.array([-3,-3]))  # 6 time TBP points
    lines[7].set_data(pdata[2], pdata[3][1])  # 7 profile compressed
    lines[8].set_data(TBPtcomp, np.array([-3,-3]))  # 8 profile compressed TBP points
    lines[11].set_data(npp, Pp / np.max(Pp))

    il = 14

    ddst = 0
    for iu in range(len(Ip)):
        lines[il + iu].set_data(np.array(Ip[iu][1]+ddst), np.array(Ip[iu][2]))
        ddst = ddst+Ip[iu][1][-1]
        #print(Ip[iu][1])
        #print(Ip[iu][2])

    nalmff1 = np.angle(yfnalm)
    nalmfb1 = np.angle(ybnalm)

    lines[33].set_data(pdata[2], np.unwrap(nalmff1 - nalmfb1)/np.pi)


    """
    # print [-2] element filterB params befor and after
    print(Ip[-2][2])
    print(Ep[-2][2])
    print(Ip[-2][0])
    print(Ep[-2][0])
    """
    """
    line41, = axs[3].plot([], [])  # 14 passive fiber energy
    line42, = axs[3].plot([], [])  # 15 coupler energy
    line43, = axs[3].plot([], [])  # 16 active fiber energy
    line44, = axs[3].plot([], [])  # 17 filter energy
    line45, = axs[3].plot([], [])  # 18 sesam energy
    """


    return lines
def tcf(x2, uf2):
    # func to get Taylor coeff of data x2 - uf2.
    # Here the set of linear eq for Taylor coeff has been got
    # R.dot(Taylor coeff) = uf2
    # Taylor coeff = inv(R).dot(uf2)
    R = np.zeros([1,len(x2)-1])
    for zi in range(1, len(x2)):
        dR = []
        for iu in range(1, len(x2)):
            d = (x2[zi] - x2[0]) ** iu / math.factorial(iu)
            dR = np.hstack([dR, d])
        R = np.vstack([R, dR])
    R = R[1:,:]

    iR = np.linalg.inv(R)
    Xt = iR.dot(uf2[1:] - uf2[0])
    return Xt
def tvl(Xt, x3, y0, x30):
    # Egz = tvl(Xt, z, Ez[0], La[0])
    dr = 0
    for zi in range(0, len(Xt)):
        dr1 = Xt[zi] * (x3 - x30) ** (zi + 1) / math.factorial((zi + 1))
        dr = dr + dr1
    dr = dr + y0
    return dr
def cavityrun(ig, y0w):

    Ip = [] # data for pulse intensity through cavity
    Ep = []  # data for pulse energy through cavity
    def lzrcrcl(y0w):
        elementtogetpulseparams = 'nalmA'
        Ip = []  # data for pulse intensity through cavity
        Ep = []  # data for pulse energy through cavity
        y = y0w
        print('g0 = ', cvtstruct[0][5])
        for m in cvtstruct:
            if m[0].__name__ == 'active':
                m[9] = 1
            print(m[1], end=" ")
            start = time.time()
            y11 = m[0](y, m[1:], ig)
            end = time.time()
            print('----time---', round(end - start, 2), 's')
            y = y11[0]
            Ip.append(y11[1])
            Ep.append(y11[2])
            if m[0].__name__ == 'active':
                m[6] = y11[2][2]*0
            if m[1] == elementtogetpulseparams:
                yout = y
                print(elementtogetpulseparams, 'after pulse params=', ['Ip=', y11[1][2][1], 'Ep=', y11[2][2][1]])
            #if m[0].__name__ == 'nalm':
            #    y = yl
            #yl = y
        Pp0 = np.max(np.abs(ifft(yout))**2)

        return [y, yout, Pp0, Ip, Ep]

    def lzrlnr(y0w):
        elementtogetpulseparams = 'couplerA'
        cvtstructb = cvtstruct[-1::-1]

        Ip = []  # data for pulse intensity through cavity
        Ep = []  # data for pulse energy through cavity
        y = y0w
        print('g0 = ', getgf(somepar[3], ig))
        for m in cvtstruct:
            print(m[1]) # print cavity element name
            if m[0].__name__ == 'active':
                m[9] = 1
            y11 = m[0](y, m[1:], ig)
            y = y11[0]
            Ip.append(y11[1])
            Ep.append(y11[2])
            if m[0].__name__ == 'active':
                m[6] = y11[2][2]
            if m[1] == elementtogetpulseparams:
                yout = y
                print(['Ip=', y11[1][2][1], 'Ep=', y11[2][2][1]])
        Pp0 = np.max(np.abs(ifft(yout))**2)

        for m in cvtstructb:
            print(m[1]) # print cavity element name
            if m[0].__name__ == 'active':
                m[9] = -1
            y11 = m[0](y, m[1:], ig)
            y = y11[0]
            Ip.append(y11[1])
            Ep.append(y11[2])
            if m[0].__name__ == 'active':
                m[6] = y11[2][2]
            #if m[1] == elementtogetpulseparams:
            #    yout = y
            #    print(['Ip=', y11[1][2][1], 'Ep=', y11[2][2][1]])
        #Pp0 = np.max(np.abs(ifft(yout))**2)

        return [y, yout, Pp0, Ip, Ep]

    cavitytype = {1: lzrcrcl, 2: lzrlnr}
    yrep, yout, Pp0, Ip, Ep = cavitytype[1](y0w)



    return yrep, yout, Pp0, Ip, Ep
def movetotime(yw,tam):
    yt = np.abs(ifft(yw))
    ytl = yt.tolist()
    po = max(ytl)
    ipo = ytl.index(po)
    res = t[ipo]-tam
    ywmoved = yw*np.exp(1j * w * res)
    return ywmoved
def compresspulse(x, yw):


    #d0 = -1 * x * w0 ^ 3 * Lg ^ 2 * cosTr ^ 3 / 4 / pi ^ 2 / c;
    #bg = 4 * pi ^ 2 * c * d0 / w0 ^ 4 / Lg ^ 2 / (cosTr) ^ 5 * (1 + sinTr * sinTi) * 3 / 2;
    ut2 = ifft(yw * np.exp(1j * x * w **2))
    cpw = 1 / np.max(np.abs(ut2))

    return cpw
def getpulsecompressed(yw, x):
    uw2 = yw * np.exp(1j * x * w **2)
    return uw2
def changecvtstruct(ij):

    # nalm passive fiber L change
    L1 = 30
    L12 = 2
    ii2 = ij - 1
    if ii2 < 0:
        ii2 = 0
    L = L1 + L12*(ii2 //(multiple2savedata))
    Lnalm2 = np.arange(0, L, 0.01)
    Lnalm2 = np.hstack([Lnalm2, L])

    # nalm g0 change
    g0nalm1 = 20
    g0nalm2 = 0
    g0nalm = g0nalm1 + g0nalm2 * (ii2 // multiple2savedata)

    # cavity g0 change
    g0a1 = 20
    if ij < 20:
        g02 = (g0a1 - 0.01 * g0a1)/20 * ij +0.01 * g0a1
    else:
        g02 = g0a1

    cvtstruct[2][4][-1] = Lnalm2 # nalm passiveD fiber length update with i - cavity round trip number
    cvtstruct[2][3][5] = g0nalm # nalm activeB fiber g0 update with i
    cvtstruct[0][5] = g02 # cavity activeA fiber g0 update with i



    return 0
def cavityshow(i, ie):

    print()
    print()
    print("Iteration = ", i)

    changecvtstruct(i)
    yrep, yout, Pp1, Ip, Ep = cavityrun(i, y0w1[0])
    Pp.append(Pp1)
    ywmoved = movetotime(yrep, 0)
    lines = report(yout, Pp, Ip, Ep)
    y0w1[0] = ywmoved

    if i > 0:
        if i//multiple2savedata == i / multiple2savedata:
            savetoxls(i / multiple2savedata, yout)

    if i == ie:
        plotgraph(yout, Pp, Ip, Ep)
        workbook.close()
        print('file closed')
    return lines
def plotgraph(ywtlmd, Pp, Ip, Ep):

    global yfnalm, ybnalm

    pdata = pulseparams(ywtlmd)  # to get pdata
    gain = afgain(somepar[1])
    gain = 10 * np.log10((np.abs(gain) ** 2) / np.max(np.abs(gain) ** 2))
    gain = fftshift(gain)
    # end gain

    # get filter
    flt = filterd(somepar[2])
    flt = 10 * np.log10((np.abs(flt) ** 2) / np.max(np.abs(flt) ** 2))
    flt = fftshift(flt)
    # end filter

    #axs[3].set_xlim([0, np.max(np.array(Ip[0][1])*1.5)])
    #axs[3].set_ylim([0, np.max(np.array(Ip[2][2])*1.5)])
    #axs[2].set_xlim([0, npp[-1]+1])
    #axs[2].set_ylim([0, np.max(Pp)])
    npp = np.arange(1, len(Pp) + 1)

    TBPlmd = pdata[7]
    TBPt = pdata[6][0]
    TBPtcomp = pdata[6][1]

    fig, ax0 = plt.subplots()
    ax0.plot(pdata[0], pdata[1])  # 0 spectrum intesity log10
    ax0.plot(pdata[0], gain)  # gain 1
    ax0.plot(pdata[0], flt)  # filter 1
    ax0.plot(TBPlmd, np.array([-3,-3]), 'o')  # TBPpoints lmd

    ax0.legend(['Spectrum', 'Gain', 'Filter', '-3dB level'], loc='lower left')
    ax0.set_xlim([lmd[0], lmd[-1]])  # axs[0] - to show spectrum
    ax0.set_ylim([-100, 0.1])
    ax0.set_title('Spectrums, dB')
    ax0.set_xlabel('Wavelength, nm')

    fig, ax1 = plt.subplots()
    ax1.plot(pdata[2], pdata[3][0])  # 4 profile intencity log10
    ax1.plot(TBPt, np.array([-3,-3]),'o')  # 6 time TBP points
    ax1.plot(pdata[2], pdata[3][1])  # 7 profile compressed
    ax1.plot(TBPtcomp, np.array([-3,-3]),'o')  # 8 profile compressed TBP points
    ax1.legend(['Intensity', '-3dB Intsty', 'Int compressed', '-3dB Int cmprsd'], loc='lower left')
    axch = ax1.twinx()
    axch.plot(pdata[4], pdata[5], 'r')  # 5 chirp
    axch.legend(['Chirp, nm'], loc='upper right')
    ax1.set_ylim([-100, 0.1])
    axch.set_ylim([np.min(pdata[5]), np.max(pdata[5])])
    ax1.set_title('Pulse intensity, dB')
    ax1.set_xlabel('Time delay, ps')

    fig, ax2 = plt.subplots()
    ax2.plot(npp, Pp)

    ax2.set_title('Peak intensity, W')
    ax2.set_xlabel('Cavity round trip')

    il = 14
    ddst = 0
    fig, ax3 = plt.subplots()
    for iu in range(len(Ip)):
        ax3.plot(np.array(Ip[iu][1]+ddst), np.array(Ip[iu][2]))
        ddst = ddst+Ip[iu][1][-1]

    ax3.set_title('Peak intensity, W')
    ax3.set_xlabel('Cavity coordinate, m')

    leg3 = []
    for ity in cvtstruct:
        leg3.append(ity[1])
    ax3.legend(leg3)
    """
    il = 14

    ddst = 0
    for iu in range(len(Ip)):
        lines[il + iu].set_data(np.array(Ip[iu][1]+ddst), np.array(Ip[iu][2]))
        ddst = ddst+Ip[iu][1][-1]
    """
    nalmff1 = np.angle(yfnalm)
    nalmfb1 = np.angle(ybnalm)

    axsnalm.plot(pdata[2], np.unwrap(nalmff1 - nalmfb1)/np.pi)
    plt.grid(True)

    plt.show()
    plt.draw()
    return 0
def savetoxls(ii, yout):
    ii2 = int(ii)
    pdata = pulseparams(yout)[6:]  # to get pdata
    #[TBPt, TBPtcomp], TBPlmd, [TBPdata, TBPdatacomp], [TBPt60, TBPtcomp60], TBPlmd60, Ip, Ep
    TBPt = pdata[0][0][1] - pdata[0][0][0]
    TBPtcomp = pdata[0][1][1] - pdata[0][1][0]
    TBPlmd = pdata[1][1] - pdata[1][0]
    TBPdata = pdata[2][0]
    TBPdatacomp = pdata[2][1]
    TBPt60 = pdata[3][0][1] - pdata[3][0][0]
    TBPtcomp60 = pdata[3][1][1] - pdata[3][1][0]
    TBPlmd60 = pdata[4][1] - pdata[4][0]
    Ip = pdata[5]
    Ep = pdata[6]

    # nalm passiveD L
    CavityParamToSave = cvtstruct[2][4][-1][-1]

    worksheet.write(ii2, 0, ii2)
    worksheet.write(ii2, 1, TBPt)
    worksheet.write(ii2, 2, TBPtcomp)
    worksheet.write(ii2, 3, TBPlmd)
    worksheet.write(ii2, 4, TBPdata)
    worksheet.write(ii2, 5, TBPdatacomp)
    worksheet.write(ii2, 6, TBPt60)
    worksheet.write(ii2, 7, TBPtcomp60)
    worksheet.write(ii2, 8, TBPlmd60)
    worksheet.write(ii2, 9, Ip)
    worksheet.write(ii2, 10, Ep/1e-9)
    worksheet.write(ii2, 11, CavityParamToSave)

    print('save ii=', ii2)

    return 0
    #######


"""
La = np.linspace(-3,10,100)
Ez = 1 - 0.1 * La - 1 * La**2 + (La-5)**3
x2 = np.linspace(La[0], La[-1], 12)
uf = interp1d(La, Ez, kind='cubic')
uf2 = uf(x2)
Xt = tcf(x2, uf2) # Taylor coeff for Ez profile through active fiber length
Egz2 = []
for z in x2:
    Egz = tvl(Xt, z, Ez[0], La[0])
    Egz2.append(Egz)
plt.plot(La, Ez)
plt.plot(x2, Egz2, 'o')
plt.show()
plt.draw()
dddddddd
"""


workbook = xlsxwriter.Workbook('Pulse.xlsx')
worksheet = workbook.add_worksheet()
initxlsfile()

multiple2savedata = 50 # cavity round trip multiple to save pulse data to file

c = 3 * 1e8 / 1e12 # speed of light in m/ps
lmc = 1030 # active fiber amplification band center
w0 = 2 * np.pi * c / lmc / 1e-9

twidth = 1000  # time domain width in ps
n1 = 2 ** 11  # number of points in time domain scale

frames0 = 4000 # number of cavity runs

y0w, t, dt, w, lmd = initdata(twidth, n1) # to get: y0w - seed pulse spectrum, t - time domain points vector, dt - time interval, w - angular frequency of considered spectrum, lmd - pulse spectrum width in wavelength
cvtstruct, somepar = cavitystruct() # cor - list of names of cavity elements, cel - dict of cavity elements parameters

disptocompress0 = [-2] # initial dispertion to compress pulse
y0w1 = [y0w] # initial pulse spectrum
Pp = [] # sequence of pulse intensity to see convergence



global yfnalm, ybnalm

# graph init
fig, axs = plt.subplots(1, 4)  # graphical elements to show results
axch = axs[1].twinx()

fignalm, axsnalm = plt.subplots()
plt.grid(True)

lines = initgraph()

# FuncAnimation launch simulation: cavityshow -> cavityrun -> ...
anim = FuncAnimation(fig, cavityshow, fargs=(frames0-1,), init_func=initgraph, frames=frames0, interval=1000, blit=True, repeat=False) # function to show axs

mp.show()
mp.draw()
"""
fig, axs = plt.subplots(1, 4)
axs[0].set_xlim([lmd[0], lmd[-1]])
axs[0].set_ylim([-100, 0.1])
axs[1].set_xlim([t[0], t[-1]])
axs[1].set_ylim([-100, 0.1])
axs[2].set_xlim([0, 1500])
axs[2].set_ylim([0, 50])
axs[3].set_xlim([0, 1.5])
axs[3].set_ylim([0, 3])
lines = cavityshow(1)

#report(yrep, cp)
#print(yrep)
mp.show()
mp.draw()

"""

"""
fig = plt.figure()
ax = plt.axes(xlim=(0, 4), ylim=(-2, 2))
line, = ax.plot([], [], lw=3)
line.set_data(t, y)
"""


