#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# Load data from CSV file
data = np.loadtxt('runtime_diags.csv', delimiter=',', skiprows=1)

# Extract columns based on CSV structure:
# Column 1: time
# Column 8: KE.O+
# Column 14: KE.e-
# Column 15: PE
# Column 16: E_total
time = data[:, 1]           # time column
KE_oxygen = data[:, 8]      # KE.O+ 
KE_electrons = data[:, 14]  # KE.e-
PE = data[:, 15]            # PE
E_total = data[:, 16]       # E_total

# Create plot
plt.figure(figsize=(10, 6))
plt.plot(time, KE_electrons, label='KE (e⁻)', linewidth=2)
plt.plot(time, KE_oxygen, label='KE (O⁺)', linewidth=2)
plt.plot(time, PE, label='PE', linewidth=2)
plt.plot(time, E_total, label='E_total', linewidth=2, linestyle='--')

plt.xlabel('Time (s)', fontsize=12)
plt.ylabel('Energy (J)', fontsize=12)
plt.title('Energy Components vs Time', fontsize=14)
plt.legend(fontsize=10)
plt.grid(True, alpha=0.3)
plt.tight_layout()

# Save and show
plt.savefig('energy_plot.png', dpi=300)
plt.show()
