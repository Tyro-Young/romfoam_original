import matplotlib

matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt

from smt.sampling_methods import LHS

xlimits = np.array([[-1, 1], [-1, 1], [-1, 1]])
sampling = LHS(xlimits=xlimits, criterion="cm")

num = 15
x = sampling(num)

print((x.shape))
print(x)
plt.plot(x[:, 0], x[:, 1], "o")
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("LHS.png", bbox_inches="tight")  # save the figure to file

np.savetxt("twistSamples.csv", x, delimiter=",")

