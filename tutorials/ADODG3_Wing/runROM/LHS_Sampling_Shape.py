import matplotlib

matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt

from smt.sampling_methods import LHS

xlimits = np.array(
    [
        [-0.01, 0.01],
        [-0.01, 0.01],
        [-0.01, 0.01],
        [-0.01, 0.01],
        [-0.01, 0.01],
        [-0.01, 0.01],
        [-0.01, 0.01],
        [-0.01, 0.01],
        [-0.01, 0.01],
        [-0.01, 0.01],
        [-0.01, 0.01],
        [-0.01, 0.01],
        [-0.01, 0.01],
        [-0.01, 0.01],
        [-0.01, 0.01],
        [-0.01, 0.01],
        [-0.01, 0.01],
        [-0.01, 0.01],
    ]
)
sampling = LHS(xlimits=xlimits, criterion="cm")

num = 90
x = sampling(num)

print((x.shape))
print(x)
plt.plot(x[:, 0], x[:, 1], "o")
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("LHS.png", bbox_inches="tight")  # save the figure to file

np.savetxt('shapeSamples.csv', x, delimiter=',')

