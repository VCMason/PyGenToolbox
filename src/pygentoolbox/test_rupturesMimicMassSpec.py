import matplotlib.pyplot as plt
import ruptures as rpt

print('Keep the number of genes < 100 for plotting. Note: time points are rows in the data matrix but columns in the plot. genes are columns in data matrix but rows in plot\n')
rows = int(input('How many time points would you like to generate in the data matrix? (ex: 15) : '))
columns = int(input('How many genes would you like to generate in the data matrix? (ex: 50) : '))
display = input('Do you want the plot? (y/n): ')
n_bkps = int(input('How many break points? (ex: 1) : '))
m = input('Which Model : l1, l2, or rbf ? : ')
p = int(input('penalty for model? ex: (integers 1-10) : '))

# generate signal
n_samples, dim, sigma = rows, columns, 4
#n_bkps = 1  # number of breakpoints
signal, bkps = rpt.pw_constant(n_samples, dim, n_bkps, noise_std=sigma)

# detection
algo = rpt.Pelt(model=m).fit(signal)
result = algo.predict(pen=p)

# display # suppress this display if you want to test more 
if display == 'y':
    rpt.display(signal, bkps, result)
    plt.show()

print('rows: %d columns: %d' % (len(signal), len(signal[0])))
print('Input data format:')
print(signal)
print('Generated Break Points')
print(bkps)

print('Result:')
print(result)