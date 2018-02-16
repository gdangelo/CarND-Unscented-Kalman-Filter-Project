import matplotlib.pyplot as plt
import numpy as np

def read_file(filepath):
    data = []
    with open(filepath) as f:
        str_data = f.read()
        data = str_data.split()
    return data

def plot_graph(x, y, threshold, sensor_type, percentage):
    # Clear figure
    plt.clf()

    fig = plt.figure()

    # NIS values
    plt.plot(x, y, color='b', label='NIS values at each step')
    # In 5% of the time nis values should be higher than threshold
    plt.axhline(y=threshold, color='r', label='threshold (' + str(threshold) + ')')

    # Add axis labels
    plt.xlabel('step number')
    plt.ylabel('nis')

    # Add titles
    title = 'Normalize Innovation Squared (NIS) graph'
    subtitle = '{:.2f}'.format(percentage) + '% NIS values above threshold'
    plt.title(title + '\n' + subtitle)

    plt.legend(loc='upper right')

    # Save graph
    if (sensor_type == 'L'):
        fig.savefig('nis_lidar_graph.png')
    elif (sensor_type == 'R'):
        fig.savefig('nis_radar_graph.png')

def main():
    # Read NIS values from text files for both Radar & LIDAR
    nis_laser = read_file("NIS_Lidar_.txt")
    nis_radar = read_file("NIS_Radar_.txt")

    # Store the values in a numpy array
    nis_laser = np.array(nis_laser, dtype='f')
    nis_radar = np.array(nis_radar, dtype='f')

    # Find percentage above threshold
    laser_threshold = 5.99
    radar_threshold = 7.81
    nis_above_laser_thres = len(np.where(nis_laser > laser_threshold)[0]) / float(len(nis_laser)) * 100
    nis_above_radar_thres = len(np.where(nis_radar > radar_threshold)[0]) / float(len(nis_radar)) * 100

    # Plot the NIS graph using matplotlib
    plot_graph(range(len(nis_laser)), nis_laser, laser_threshold, 'L', nis_above_laser_thres)
    plot_graph(range(len(nis_radar)), nis_radar, radar_threshold, 'R', nis_above_radar_thres)

    print("Done.")

if __name__ == "__main__":
    main();
