Momentum Algorithm

Hourly and daily surface water temperature simulations were performed with the momentum algorithm. Momentum is an iterative first order optimization method that uses gradient calculated from the average loss of a neural network. Momentum accumulates information about the gradients during the training in the way that helps to gradually forget information about previous gradients over time.

wn =wn-1 +En
wn = wn-1-wn

Where wnis a vector that stores parameters of the neural network after n training iterations, wnis a weight increment on the n-th training iteration, is a learning rate, is a parameter that controls the rate at which information about previous gradients is being erased and En is a gradient calculated during the n-th epoch.

In general, momentum converges much faster than the gradient descent algorithm. It happens because momentum manages to adopt the learning rate during the training per each parameter individually. Intuitively, when gradient vector doesn’t change its direction during the training than magnitude of the E will increase. This increase is equivalent to increase in the learning rate when E is being kept unchanged.
