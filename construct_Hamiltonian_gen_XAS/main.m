%This block generates the ligand field potential Hamiltonian in many-body
%basis
%Tight-binding Hamiltonian is constructed manually
%Index ordering is found in 'qn', defined in the next block
fprintf('Constructing tight-binding Hamiltonian')
load('H_TB_NiO')
H_TB = H_TB_NiO;
numorbitals = 52;
numelectrons = 50;
%Convert 'H_TB' in single-particle to many-body basis 'HTotal_TB'
%'binaryindices' is list of orbital occupations for each row/column of
%'HTotal_TB'
[HTotal_TB,binaryindices]=convert_HTBtoHTotal(H_TB,...
    numorbitals,numelectrons);
fprintf('Completed tight-binding Hamiltonian')

%This block generates the spin-orbit coupling Hamiltonian in many-body
%basis
fprintf('Constructing spin-orbit coupling Hamiltonian')
lambda_3d = 0.08; %spin-orbit coupling constants for Ni 3d in eV
lambda_2p = 11.51; %spin-orbit coupling constants for Ni 2p in eV
%'qn' is a list of quantum numbers where the row index corresponds
%to the same row/column index as all of the Hamiltonians in single-particle
%basis. The five columns correspond to, in order, n, l, m, spin (1 for up,
%0 for down), and atom type (1 for metal center, 0 for ligand)
load('qn');
H_SO = gen_HSO(numorbitals, qn, lambda_3d, lambda_2p);
HTotal_SO = convert_HSOtoHTotal(H_SO, numorbitals, numelectrons);
fprintf('Completed spin-orbit coupling Hamiltonian')

%This block generates the exact Coulomb correlation Hamiltonian in
%many-body basis
%'Rk_NiO' contains Slater integral values. First four columns are nlms
%quantum numbers of R1, second four columns are nlms quantum numbers of R2,
%9th column is atom type of R1 (1 for metal, 0 for ligand), 10th column is
%atom type of R2, 11th column is the multipole index k, and 12th column is
%Rk itself
%Rk values are from angular integration conducted in folder
%'gen_radialfunctions'
fprintf('Constructing exact Coulomb Hamiltonian')
load('Rk_NiO')
Rk = Rk_NiO;
%'gauntcoefficients' first four columns: l1, m1, l2, m2. Fifth column is
%atom type (1 for metal, 0 for ligand), and fifth column is Gaunt
%coefficient itself taken from Slater's textbook 'Quantum Theory of Solids'
load('gauntcoefficients')
%'V' contains value of V(alpha,beta,gamma,delta) in fifth column
V = gen_V(numorbitals, gauntcoefficients, qn, Rk);
HTotal_Coulomb = convert_HCoulombtoHTotal(V, numorbitals, numelectrons);
fprintf('Completed exact Coulomb Hamiltonian')

%This block generates the mean-field Coulomb correlation Hamiltonian in
%many-body basis
%'rho' contains expectation value of occupation of each spin orbital
%Index ordering can be found in 'qn'
fprintf('Constructing mean-field Hamiltonian')
rho = gen_rho(H_TB, numelectrons);
HSI = gen_HSelfInteraction(numorbitals, numelectrons, V, rho);
HTotal_MF = convert_HCoulombMeanFieldtoHTotal(numorbitals, numelectrons,...
    V, rho)+ convert_HSItoHTotal(HSI,numorbitals,numelectrons);
fprintf('Completed mean-field Hamiltonian')

%Construct total Hamiltonian
HTotal = HTotal_TB + HTotal_Coulomb + HTotal_SO - HTotal_MF;

%Generate XAS spectra where w is energy in eV and K is intensity
fprintf('Generating XAS spectra')
numinitialstates = 1;
nvector = [1 1 1]; %x, y, z polarization of incoming radiation
FWHM = 0.1; %full width half max of Lorentzian broadening in eV
[w, K] = gen_XASspectra(HTotal, numinitialstates,nvector,qn,...
    gauntcoefficients,FWHM)
plot(w,K)
fprintf('Completed XAS spectra')
