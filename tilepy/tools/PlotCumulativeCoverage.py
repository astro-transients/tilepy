import pandas as pd
import glob
import numpy as np
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

def ObtainCumulativeProbabilityPlot(folder_path, event_name, WhatToPlot):

    # Get a list of all .txt files in the specified folder
    
    #file_list = glob.glob(folder_path)
    #print(folder_path)

    # Initialize an empty DataFrame to store the cumulative 5th columns
    file_path =  folder_path 
    cumulative_column = None

    # Initialize a list to store individual 5th columns for plotting
    individual_columns = []

    # Loop through each file and add the 5th column to the cumulative DataFrame
    # Loop through each file and add the 5th column to the cumulative DataFrame
    # Read the .txt file into a Pandas DataFrame
    df = pd.read_csv(file_path, sep=' ')  # Adjust the separator if needed

    observatory_column = df.columns[-1]
    
    grouped_data = df.groupby(observatory_column)

    # Find the row with the maximum value in 'value_column' for each group in 'group_column'
    df['cumsum'] = df[WhatToPlot].cumsum()  # Cumulative sum of PGal values

    max_value_indices = df.groupby("Observation Time UTC")['cumsum'].idxmax()

    # Select the rows with the maximum values
    max_rows = df.loc[max_value_indices]

    # Display the resulting DataFrame
    print(max_rows)
    print('---------------')
    print('---------------')
    print(df[WhatToPlot])
    print('---------------')
    print('---------------')
    #x_values2 = pd.to_datetime(max_rows["Observation Time UTC"])
    max_rows['date_column'] = pd.to_datetime(max_rows["Observation Time UTC"])
    
    y_values2 = max_rows["cumsum"]
    # Plot the cumulative values for all observatories
    plt.plot(max_rows['date_column'], y_values2, label='All Observatories', color='black')
    

    for observatory, observatory_data in grouped_data:
        #print(f"Observatory: {observatory}")
        #print(observatory_data)
        print("\n")
        # Extract relevant columns for the plot
        x_values = pd.to_datetime(observatory_data["Observation Time UTC"])
        #print(x_values)
        y_values = observatory_data[WhatToPlot].cumsum()

        # Plot the values for each observatory
        plt.plot(x_values, y_values,label=observatory)

    # Add labels and legend
    # Format the x-axis to display date and time
    # Set the locator and formatter for the x-axis to show dates
    # Manually set date ticks on the x-axis
    # Format the x-axis to display date and time
    plt.gca().xaxis.set_major_locator(mdates.AutoDateLocator())
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    plt.gcf().autofmt_xdate()  # Auto-format the date labels for better readability

    #plt.tight_layout()  # Adjust layout for better spacing
    plt.title("%s trigger coverage on 2023-10-12" %(event_name))
    plt.xlabel('Date [UTC]')
    plt.ylabel('Total Covered Probability [%]')
    #plt.xscale('log')  # Set X-axis to log scale
    #plt.yscale('log')  # Set Y-axis to log scale
    plt.xticks(rotation=45)  # Rotate x-axis labels for better readability

    plt.legend()
    plt.grid()
    plt.savefig('%s_%s_cumulativeDistribution.png' %(event_name,WhatToPlot), dpi=300, bbox_inches='tight')
    plt.close()
    #plt.show()

def ObtainCumulativeProbabilityPlotLog(folder_path, event_name, WhatToPlot,dateToSubtract):

    # Get a list of all .txt files in the specified folder
    
    #file_list = glob.glob(folder_path)
    #print(folder_path)

    # Initialize an empty DataFrame to store the cumulative 5th columns
    file_path =  folder_path 
    cumulative_column = None

    # Initialize a list to store individual 5th columns for plotting
    individual_columns = []

    # Loop through each file and add the 5th column to the cumulative DataFrame
    # Loop through each file and add the 5th column to the cumulative DataFrame
    # Read the .txt file into a Pandas DataFrame
    df = pd.read_csv(file_path, sep=' ')  # Adjust the separator if needed

    observatory_column = df.columns[-1]
    
    grouped_data = df.groupby(observatory_column)

    # Find the row with the maximum value in 'value_column' for each group in 'group_column'
    df['cumsum'] = df[WhatToPlot].cumsum()  # Cumulative sum of PGal values

    max_value_indices = df.groupby("Observation Time UTC")['cumsum'].idxmax()

    # Select the rows with the maximum values
    max_rows = df.loc[max_value_indices]

    # Display the resulting DataFrame
    print(max_rows)
    print('---------------')
    print('---------------')
    print(df[WhatToPlot])
    print('---------------')
    print('---------------')
    dateToSubtract = pd.to_datetime(dateToSubtract)
    #x_values2 = pd.to_datetime(max_rows["Observation Time UTC"])
    max_rows['date_column'] = pd.to_datetime(max_rows["Observation Time UTC"])
    max_rows['difference'] = max_rows['date_column'] - dateToSubtract
    # Convert timedelta values to seconds
    max_rows['difference_seconds'] = max_rows['difference'].dt.total_seconds()
    y_values2 = max_rows["cumsum"]
    # Plot the cumulative values for all observatories
    print(max_rows['difference_seconds'])
    plt.plot(max_rows['difference_seconds'], y_values2, label='All Observatories', color='black')

    

    for observatory, observatory_data in grouped_data:
        #print(f"Observatory: {observatory}")
        #print(observatory_data)
        print("\n")
        # Extract relevant columns for the plot

        observatory_data['date_column'] = pd.to_datetime(observatory_data["Observation Time UTC"])
        observatory_data['difference'] = observatory_data['date_column'] - dateToSubtract
        # Convert timedelta values to seconds
        observatory_data['difference_seconds'] = observatory_data['difference'].dt.total_seconds()
        #y_values = observatory_data["cumsum"]
        y_values = observatory_data[WhatToPlot].cumsum()

        # Plot the values for each observatory

        plt.plot(observatory_data['difference_seconds'], y_values,label=observatory)

    # Add labels and legend
    # Format the x-axis to display date and time
    # Set the locator and formatter for the x-axis to show dates
    # Manually set date ticks on the x-axis
    # Format the x-axis to display date and time
    #plt.gca().xaxis.set_major_locator(mdates.AutoDateLocator())
    #plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    #plt.gcf().autofmt_xdate()  # Auto-format the date labels for better readability

    #plt.tight_layout()  # Adjust layout for better spacing
    plt.title("%s trigger coverage on 2023-10-12" %(event_name))
    plt.xlabel('T-T0 [seconds]')
    plt.ylabel('Total Covered Probability [%]')
    plt.xscale('log')  # Set X-axis to log scale
    #plt.yscale('log')  # Set Y-axis to log scale
    plt.xticks(rotation=45)  # Rotate x-axis labels for better readability

    plt.legend()
    plt.grid()
    plt.savefig('%s_%s_cumulativeDistribution.png' %(event_name,WhatToPlot), dpi=300, bbox_inches='tight')
    plt.close()
    #plt.show()