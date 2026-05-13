import numpy as np


def predict(x, y, x_predict):
    n_points = len(x)

    print("extrapolating to", x_predict, "from", x)

    numerator_matrix = np.tile(x_predict - x, (n_points, 1))
    np.fill_diagonal(numerator_matrix, 1.0)

    denominator_matrix = x[:, None] - x[None, :]
    np.fill_diagonal(denominator_matrix, 1.0)

    num_terms = np.prod(numerator_matrix, axis=1)
    denom_terms = np.prod(denominator_matrix, axis=1)

    return np.sum(num_terms * y / denom_terms)

def curvature(x, y):
    dx = x[1:] - x[:-1]
    dydx = (y[1:] - y[:-1]) / dx

    dx_mean = 0.5*(dx[1:] + dx[:-1])
    d2ydx2 = (dydx[1:] - dydx[:-1]) / dx_mean

    dydx = np.concatenate(([dydx[0]], 0.5*(dydx[1:] + dydx[:-1]), [dydx[-1]]))
    d2ydx2 = np.concatenate(([d2ydx2[0]], d2ydx2, [d2ydx2[-1]]))

    k = d2ydx2 / ((1+dydx**2)**(3/2))

    return k


def continuity_score(x, y):
    """ Scores the degree of continuity across gaps """
    pass

def sorted_data(*args):
    data = np.array(args)
    data = np.sort(data, axis=0)
    return data


if __name__ == "__main__":

    test_predictions = False
    test_curvature = True

    if test_curvature:
        n_points = 101

        x = np.linspace(0, 2 * np.pi, n_points)
        y = [1+np.sin(x), 1+np.sin(x+0.3), 1+np.sin(2*x), 0.8+np.cos(x)]
        y = y + [-yy for yy in y]

        import matplotlib.pyplot as plt

        plt.figure("Input curves")

        for yy in y:
            plt.plot(x, yy)

        plt.figure("Sorted curves")
        y = sorted_data(*y)

        for i in range(y.shape[0]):
            plt.plot(x, y[i, :])

        plt.figure("Curvature")


        for i in range(y.shape[0]):

            k = curvature(x, y[i, :])

            plt.plot(x, k)



        plt.show()



    if test_predictions:

        n_points = 101

        x = np.linspace(0, 2*np.pi, n_points)
        y = np.sin(x)
        z = np.sin(x+0.3)

        n = 4
        predictions = []
        for i in range(n, n_points):
            y_pred = predict(x[i-n:i], y[i-n:i], x[i])

            predictions.append(y_pred)

        import matplotlib.pyplot as plt

        plt.plot(x, y)
        plt.scatter(x[n:], predictions)

        plt.show()

