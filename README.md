# flavoured-cptk

Added new sections in sh_main.F and sh_init.F to implement a static oulomb barrier and XT-CT offset

Two new keywords in CT_BLOCK of input file: 

- NUMBER_SITES = total number of CT states in the system

- XT_CT_OFFSET = energy offset between interfacial CT and XT states (in Hartrees)

- include FRZ_DIAGONALS.include file below the NUMBER_SITES keyword, FRZ diagonals has the coulomb barrier value of each CT state alongside the corresponding indices of each site
