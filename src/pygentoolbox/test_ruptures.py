import matplotlib.pyplot as plt
import ruptures as rpt

# generate signal
n_samples, dim, sigma = 1000, 3, 4
n_bkps = 4  # number of breakpoints
signal, bkps = rpt.pw_constant(n_samples, dim, n_bkps, noise_std=sigma)

# detection
algo = rpt.Pelt(model="rbf").fit(signal)
result = algo.predict(pen=10)

# display
rpt.display(signal, bkps, result)
plt.show()

print('rows: %d columns: %d' % (len(signal), len(signal[0])))
print('Input data format:')
print(type(signal))
print(signal)
print('Generated Break Points')
print(bkps)

print('Result:')
print(result)