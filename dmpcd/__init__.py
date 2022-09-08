#import pcd
#import baseline
import csv
import matplotlib.pyplot as plt


def visualize_results_2d(vert_filename, edge_filename):
    x = []
    y = []
    # print('reading coords')
    with open(vert_filename, 'r') as vert_file:
        reader = csv.reader(vert_file, delimiter=' ')
        for row in reader:
            x.append(float(row[0]))
            y.append(float(row[1]))
        vert_file.close()

    edges = []
    edge_freq = {}
    morse_verts = set()
    # print('reading morse')
    with open(edge_filename) as edge_file:
        reader = csv.reader(edge_file, delimiter=' ')
        for row in reader:
            e = (int(row[0]), int(row[1]), -1)
            if e not in edge_freq:
                edge_freq[e] = 1
            else:
                # print('read all about it')
                edge_freq[e] += 1
            edges.append(e)
            # edges.append((int(row[0]), int(row[1]), int(row[2])))
            morse_verts.add(int(row[0]))
            morse_verts.add(int(row[1]))
        edge_file.close()

    morse_x = [x[i] for i in range(len(x)) if i in morse_verts]
    morse_y = [y[i] for i in range(len(y)) if i in morse_verts]
    nm_x = [x[i] for i in range(len(x)) if i not in morse_verts]
    nm_y = [y[i] for i in range(len(y)) if i not in morse_verts]

    for e in edges:
        v0 = min(e[0], e[1])
        v1 = max(e[0], e[1])
        x_vals = [x[v0], x[v1]]
        y_vals = [y[v0], y[v1]]

        plt.plot(x_vals, y_vals, c='red')

    plt.scatter(morse_x, morse_y, s=5, c='black')
    plt.scatter(nm_x, nm_y, s=5, c='black')

    plt.show()
