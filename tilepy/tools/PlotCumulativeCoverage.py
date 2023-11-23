import pandas as pd
import glob
import numpy as np
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt

def ObtainCumulativeProbabilityPlot(folder_path)

    # Get a list of all .txt files in the specified folder
    file_list = glob.glob(folder_path)

    # Initialize an empty DataFrame to store the cumulative 5th columns
    cumulative_column = None

    # Initialize a list to store individual 5th columns for plotting
    individual_columns = []

    # Loop through each file and add the 5th column to the cumulative DataFrame
    for file_path in file_list:
        # Read the .txt file into a Pandas DataFrame
        df = pd.read_csv(file_path, sep=' ', skiprows=1)  # Adjust the separator if needed

        # Extract the 5th column
        fifth_column = df.iloc[:, 3]*100  # Assumes the 5th column is at index 4 (0-based)
        cum_fifth_column = np.cumsum(fifth_column)
        # Add the 5th column to the cumulative DataFrame
        if cumulative_column is None:
            cumulative_column = cum_fifth_column
        else:
            cumulative_column = cumulative_column+cum_fifth_column

        # Store individual 5th columns for plotting
        individual_columns.append(cum_fifth_column)

    # Plot individual 5th columns
    plt.figure(figsize=(7, 6))
    for i, column in enumerate(individual_columns, start=1):
        plt.plot(column, label=f'LST-{i}')

    # Plot the cumulative 5th column
    plt.plot(cumulative_column, label='Total LST1-4', linestyle='--', linewidth=2)

    # Add labels and legend
    plt.xlabel('Observation number [#]')
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.ylabel('Total Covered Probability [%]')
    plt.legend()
    plt.grid()
    # Show the plot
    #plt.show()
    plt.savefig(folder_path+'/cumulativeDistribution.png')
