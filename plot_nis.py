import matplotlib.pyplot as plt
import numpy as np

def read_file(filepath):
    data = []
    with open(filepath) as f:
        str_data = f.read()
        data = str_data.split()
    return data

def plot_graph(x, y, threshold, sensor_type):
    plt.clf() # clear figure

    fig = plt.figure()

    plt.xlabel('step number')
    plt.ylabel('nis')
    plt.title('Normalize Innovation Squared (NIS) graph')

    plt.plot(x, y, color='b') # nis values
    plt.axhline(y=threshold, color='r') # in 5% of the time nis values will be higher than threshold

    if (sensor_type == 'L'):
        fig.savefig('nis_lidar_graph.png')
    elif (sensor_type == 'R'):
        fig.savefig('nis_radar_graph.png')

def main():
    # Read NIS values from text files for both Radar & LIDAR
    nis_laser = read_file("NIS_Lidar_.txt")
    nis_radar = read_file("NIS_Radar_.txt")

    # Store the values in a numpy array
    nis_laser = np.array(nis_laser)
    nis_radar = np.array(nis_radar)

    # Plot the NIS graph using matplotlib
    plot_graph(range(len(nis_laser)), nis_laser, 5.99, 'L')
    plot_graph(range(len(nis_radar)), nis_radar, 7.81, 'R')

    print("Done.")

if __name__ == "__main__":
    main();
