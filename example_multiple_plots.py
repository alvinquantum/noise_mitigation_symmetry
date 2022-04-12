import matplotlib.pyplot as plt

f, ax = plt.subplots(2,1, figsize=(6,6), sharex=True)

ax[0].plot([1,2,3])
ax[1].plot([3,2,1], label='Hi')
ax[1].legend()
plt.subplots_adjust(hspace=0.05) 
