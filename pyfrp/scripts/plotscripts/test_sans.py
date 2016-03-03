from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Bitstream Vera Sans', 'Lucida Grande', 'Verdana', 'Geneva']
#, Lucid, Arial, Helvetica, Avant Garde, sans-serif']
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
ax.plot([1,2,3], label='test')

ax.legend()
plt.show()