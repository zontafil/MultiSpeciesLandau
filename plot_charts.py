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
    def plotTemperatureSpeciesAxes(self, TspeciesAxes, enableLogScale, enableAnalytic):
        self.plt.clf()
        times = self.times.copy()
        self.plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

        for s in range(0, self.nspecies):
            self.plt.scatter(times, TspeciesAxes[s, 0], label=self.specieNames[s]+" x", s=0.5)
            self.plt.plot(times, TspeciesAxes[s, 0])
            self.plt.scatter(times, TspeciesAxes[s, 1], label=self.specieNames[s]+" y", s=0.5)
            self.plt.plot(times, TspeciesAxes[s, 1])

        # analytic solution
        if enableAnalytic:
            Tanalytic = Tisothermal(times, self.ma, self.mb, self.na, self.Tax, self.Tay, self.Tbx, self.Tby, self.clog00, self.clog01, self.clog11)
            self.plt.scatter(times, Tanalytic[0], color="black", label="Analytic", s=0.5)
            self.plt.plot(times, Tanalytic[0], color="black", linestyle="dashed")
            self.plt.scatter(times, Tanalytic[1], s=0.5)
            self.plt.plot(times, Tanalytic[1], color="black", linestyle="dashed")
            self.plt.scatter(times, Tanalytic[2], s=0.5)
            self.plt.plot(times, Tanalytic[2], color="black", linestyle="dashed")
            self.plt.scatter(times, Tanalytic[3], s=0.5)
            self.plt.plot(times, Tanalytic[3], color="black", linestyle="dashed")
        legend = self.plt.legend(loc='upper right')
        legend = self.plt.legend(bbox_to_anchor=(1,1.15), loc="upper right")

        for s in range(0, 2*self.nspecies):
            legend.legendHandles[s]._sizes = [30]

        if enableLogScale:
            self.plt.xscale("log")

        self.plt.grid()
        self.plt.xlabel("Time [s]")
        self.plt.ylabel("Temperature [eV]")
        self.plt.gcf().subplots_adjust(left=0.15)
        current_values = self.plt.gca().get_yticks()
        self.plt.gca().set_yticklabels(['{:.0f}'.format(x) for x in current_values])
        self.plt.savefig("out/TemperatureSpeciesAxes.eps")
        self.plt.savefig("out/TemperatureSpeciesAxes.png")

    def plotTemperatureSpecies(self, Tspecies, enableLogScale, enableAnalytic, Te1):
        # plot single species temeperature
        self.plt.clf()
        self.plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        for s in range(0, self.nspecies):
            self.plt.scatter(self.times, Tspecies[s], label=self.specieNames[s], s=0.5)
            self.plt.plot(self.times, Tspecies[s])

        if enableAnalytic:
            self.plt.axhline(y=Te1, color="black", linestyle="dashed", label="Analytic")

        legend = self.plt.legend()
        for s in range(0, self.nspecies):
            legend.legendHandles[s]._sizes = [30]
        self.plt.xlabel("Time [s]")
        self.plt.ylabel("Temperature [eV]")
        self.plt.gcf().subplots_adjust(left=0.17)

        if enableLogScale:
            self.plt.xscale("log")
        # self.plt.ylim(bottom=1130, top=1134)
        # self.plt.xlim(left=10**(-5))
        current_values = self.plt.gca().get_yticks()
        self.plt.grid()
        self.plt.savefig("out/TemperatureSpecies.eps")
        self.plt.savefig("out/TemperatureSpecies.png")

    def plotEnergy(self, Eerr):
        # plot energy error
        self.plt.clf()
        self.plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        self.plt.scatter(self.times, Eerr, s=0.5)
        self.plt.plot(self.times, Eerr)
        self.plt.xlabel("Time [s]")
        self.plt.grid()
        self.plt.ylabel("(K-K0)/K0")
        self.plt.gcf().subplots_adjust(left=0.15)
        self.plt.savefig("out/EnergyError.eps")
        self.plt.savefig("out/EnergyError.png")

    def plotEnergySpecies(self, Especies):
        # plot single species energy
        self.plt.clf()
        self.plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        for s in range(0, self.nspecies):
            self.plt.scatter(self.times, Especies[s], label=self.specieNames[s], s=0.5)
            self.plt.plot(self.times, Especies[s])
        legend = self.plt.legend()
        for s in range(0, self.nspecies):
            legend.legendHandles[s]._sizes = [30]
        self.plt.xlabel("Time [s]")
        self.plt.ylabel("Total Energy [J]")
        self.plt.grid()
        self.plt.gcf().subplots_adjust(left=0.15)
        self.plt.ticklabel_format(axis='y', style='sci', scilimits=(3,4))
        self.plt.savefig("out/EnergySpecies.eps")
        self.plt.savefig("out/EnergySpecies.png")

    def plotPnorm(self, Pnorm):
        self.plt.clf()
        self.plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        self.plt.scatter(self.times, Pnorm, s=0.5)
        self.plt.plot(self.times, Pnorm)
        self.plt.xlabel("Time [s]")
        self.plt.grid()
        self.plt.ylabel("|P|")
        self.plt.savefig("out/P.eps")
        self.plt.savefig("out/P.png")

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

    def plotVAxes(self, VFlowSpecies, enableLogScale, enableAnalytic):
        times = self.times[1:]
        flow = VFlowSpecies[:, :, 1:]
        self.plt.clf()

        self.plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        self.plt.ticklabel_format(style='sci', axis='y', scilimits=(0,2))
        for s in range(0, self.nspecies):
            self.plt.scatter(times, flow[s, 0], label=self.specieNames[s]+" x", s=0.5)
            self.plt.plot(times, flow[s, 0])
            self.plt.scatter(times, flow[s, 1], label=self.specieNames[s]+" y", s=0.5)
            self.plt.plot(times, flow[s, 1])

        if enableAnalytic:
            an = PslowDown(times, VFlowSpecies[0, 0, 0], self.ma, self.mb, self.Ta, self.Tb, self.clog00, self.na)
            self.plt.scatter(times, an, label="Analytic", s=0.5)
            self.plt.plot(times, an)

        legend = self.plt.legend()
        for s in range(0, 2*self.nspecies):
            legend.legendHandles[s]._sizes = [30]
        legend.legendHandles[2*self.nspecies]._sizes = [30]

        if enableLogScale:
            self.plt.xscale("log")

        self.plt.xlabel("Time [s]")
        self.plt.ylabel("V [m/s]")
        self.plt.grid()
        self.plt.gcf().subplots_adjust(left=0.15)
        self.plt.savefig("out/Vaxes_xlog.eps")
        self.plt.savefig("out/Vaxes_xlog.png")
