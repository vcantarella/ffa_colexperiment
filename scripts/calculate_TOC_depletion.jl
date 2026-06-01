
using Statistics

ϕ = 0.2
ρₛ = 2650 #kg/m3
L = 0.08 #column length [m]
D = 0.035 #column diameter [m]
Vt = π*D^2/4 * L #total column volume in m3
Vw = Vt*ϕ # volume of water (Pore-Volume)
Q = 2/1e6*24 #ml/hr -> m3/day
t = 30 # days
PVs = Q*t/Vw # total pore-volumes flushed
ΔC_PV = 5/4*Vw # mols of C consumed: 5/4 mmol/L in one Pore Volume
ΔC_PV*PVs # mols of C consumed during the experiment
# convert to TOC
TOC_consumed = ΔC_PV*PVs/((Vt-Vw)*ρₛ)