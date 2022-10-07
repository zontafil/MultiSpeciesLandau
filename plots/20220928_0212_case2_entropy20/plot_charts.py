from plot_analytic_times import PslowDown, Tisothermal
import numpy as np

class MultiSpeciesCharts:
    def __init__(self, plt, times, nspecies, ma, mb, Tax, Tay, Tbx, Tby, na, specieNames, clog00, clog01, clog11):
        self.plt = plt
        self.clog00 = clog00
        self.clog01 = clog01
        self.clog11 = clog11
        self.times = times
        self.nspecies = nspecies
        self.ma = ma
        self.mb = mb
        self.Tax = Tax
        self.Tay = Tay
        self.Tbx = Tbx
        self.Tby = Tby
        self.na = na
        self.specieNames = specieNames

        self.Ta = 0.5*(self.Tax + self.Tay)
        self.Tb = 0.5*(self.Tbx + self.Tby)

        # plot single species temeperature separated for each axis plt.clf()
    def plotTemperatureAnalytic(self, TspeciesAxes):
        self.plt.clf()
        times_log = self.times.copy()
        # times_log[0] = 10**(-7)
        self.plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        Tanalytic = Tisothermal(times_log, self.ma, self.mb, self.na, self.Tax, self.Tay, self.Tbx, self.Tby, self.clog00, self.clog01, self.clog11)
        for s in range(0, self.nspecies):
            self.plt.scatter(times_log, TspeciesAxes[s, 0], label=self.specieNames[s]+" x", s=0.5)
            self.plt.plot(times_log, TspeciesAxes[s, 0])
            self.plt.scatter(times_log, TspeciesAxes[s, 1], label=self.specieNames[s]+" y", s=0.5)
            self.plt.plot(times_log, TspeciesAxes[s, 1])
        # analytic solution for electrons
        self.plt.scatter(times_log, Tanalytic[0], color="black", label="Analytic", s=0.5)
        self.plt.plot(times_log, Tanalytic[0], color="black", linestyle="dashed")
        self.plt.scatter(times_log, Tanalytic[1], s=0.5)
        self.plt.plot(times_log, Tanalytic[1], color="black", linestyle="dashed")
        self.plt.scatter(times_log, Tanalytic[2], s=0.5)
        self.plt.plot(times_log, Tanalytic[2], color="black", linestyle="dashed")
        self.plt.scatter(times_log, Tanalytic[3], s=0.5)
        self.plt.plot(times_log, Tanalytic[3], color="black", linestyle="dashed")
        legend = self.plt.legend(loc='upper right')
        legend = self.plt.legend(bbox_to_anchor=(1,1.15), loc="upper right")

        for s in range(0, 5):
            legend.legendHandles[s]._sizes = [30]
        # self.plt.ylim(bottom=0)
        # self.plt.ylim(180,420)
        self.plt.ylim(250,330)
        # self.plt.ylim(250,330)
        self.plt.xlim(left=1*10**(-3), right=7*10**(-3))
        self.plt.grid()
        # self.plt.xscale("log")
        self.plt.xlabel("Time [s]")
        self.plt.ylabel("Temperature [eV]")
        self.plt.gcf().subplots_adjust(left=0.15)
        current_values = self.plt.gca().get_yticks()
        self.plt.gca().set_yticklabels(['{:.0f}'.format(x) for x in current_values])
        self.plt.savefig("out/TemperatureSpeciesAxesAnalytic_log.eps")
        self.plt.savefig("out/TemperatureSpeciesAxesAnalytic_log.png")

    def plotEnergy(self, Eerr):
        # plot energy error
        self.plt.clf()
        self.plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        self.plt.scatter(self.times, Eerr, s=0.5)
        self.plt.plot(self.times, Eerr)
        self.plt.xlabel("Time [s]")
        self.plt.grid()
        self.plt.ylabel("(E-E0)/E0")
        self.plt.gcf().subplots_adjust(left=0.15)
        self.plt.savefig("out/EnergyError.eps")
        self.plt.savefig("out/EnergyError.png")

    # plot momentum error
    def plotPerr(self, Perr, Pnorm):
        PerrNorm = np.sqrt(Perr[:,0]**2 + Perr[:,1]**2)
        PerrNorm = PerrNorm / Pnorm
        self.plt.clf()
        self.plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        self.plt.scatter(self.times, PerrNorm, s=0.5)
        self.plt.plot(self.times, PerrNorm)
        self.plt.xlabel("Time [s]")
        self.plt.ylabel("|P - P0|/|P0|")
        self.plt.grid()
        self.plt.gcf().subplots_adjust(left=0.15)
        self.plt.savefig("out/PerrNorm.eps")
        self.plt.savefig("out/PerrNorm.png")

    def plotEntropy(self, dS):
        # plot system momentum and energy
        offset = 1
        times = self.times[offset:]
        dS1 = dS[offset:]
        self.plt.clf()
        self.plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        self.plt.scatter(times, dS1, s=0.5)
        self.plt.plot(times, dS1)
        self.plt.xlabel("Time [s]")
        self.plt.grid()
        self.plt.ylabel("(S-S0)/S0")
        self.plt.xscale("log")
        # self.plt.yscale("log")
        # self.plt.ylim(bottom=0.00875, top=0.0088)
        # self.plt.ylim(bottom=0.008, top=0.0088)
        self.plt.gcf().subplots_adjust(left=0.20)
        self.plt.savefig("out/dS.eps")
        self.plt.savefig("out/dS.png")

    def plotVAxes(self, VFlowSpecies):
        self.plt.clf()
        an = PslowDown(self.times, VFlowSpecies[0, 0, 0], self.ma, self.mb, self.Ta, self.Tb, self.clog00, self.na)
        self.plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        self.plt.ticklabel_format(style='sci', axis='y', scilimits=(0,2))
        for s in range(0, self.nspecies):
            self.plt.scatter(self.times, VFlowSpecies[s, 0], label=self.specieNames[s]+" x", s=0.5)
            self.plt.plot(self.times, VFlowSpecies[s, 0])
            self.plt.scatter(self.times, VFlowSpecies[s, 1], label=self.specieNames[s]+" y", s=0.5)
            self.plt.plot(self.times, VFlowSpecies[s, 1])
        self.plt.scatter(self.times, an, label="Analytic", s=0.5)
        self.plt.plot(self.times, an)
        legend = self.plt.legend()
        for s in range(0, 2*self.nspecies):
            legend.legendHandles[s]._sizes = [30]
        legend.legendHandles[2*self.nspecies]._sizes = [30]
        self.plt.xlabel("Time [s]")
        self.plt.ylabel("V [m/s]")
        self.plt.grid()
        self.plt.gcf().subplots_adjust(left=0.15)
        self.plt.savefig("out/Vaxes.eps")
        self.plt.savefig("out/Vaxes.png")

    def plotVAxesLog(self, VFlowSpecies):
        times = self.times[1:]
        flow = VFlowSpecies[:, :, 1:]
        self.plt.clf()
        # times1 = times[1:]
        an = PslowDown(times, VFlowSpecies[0, 0, 0], self.ma, self.mb, self.Ta, self.Tb, self.clog00, self.na)
        self.plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        self.plt.ticklabel_format(style='sci', axis='y', scilimits=(0,2))
        for s in range(0, self.nspecies):
            self.plt.scatter(times, flow[s, 0], label=self.specieNames[s]+" x", s=0.5)
            self.plt.plot(times, flow[s, 0])
            self.plt.scatter(times, flow[s, 1], label=self.specieNames[s]+" y", s=0.5)
            self.plt.plot(times, flow[s, 1])
        self.plt.scatter(times, an, label="Analytic", s=0.5)
        self.plt.plot(times, an)
        legend = self.plt.legend()
        for s in range(0, 2*self.nspecies):
            legend.legendHandles[s]._sizes = [30]
        legend.legendHandles[2*self.nspecies]._sizes = [30]
        # self.plt.ylim(bottom=0)
        # self.plt.yscale("log")
        self.plt.xscale("log")
        self.plt.xlabel("Time [s]")
        self.plt.ylabel("V [m/s]")
        # self.plt.ylim([-100000,100000])
        self.plt.grid()
        self.plt.gcf().subplots_adjust(left=0.15)
        current_values = self.plt.gca().get_yticks()
        # self.plt.gca().set_yticklabels(['{:.0f}'.format(x) for x in current_values])
        self.plt.savefig("out/Vaxes_xlog.eps")
        self.plt.savefig("out/Vaxes_xlog.png")
