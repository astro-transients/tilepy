import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd

colors = ["b", "g", "r", "c", "m", "y", "k", "#1f77b4", "#ff7f0e", "#2ca02c"]


def ObtainCumulativeProbabilityPlot(
    folder_path, event_name, WhatToPlot, interactive=False
):

    # Initialize an empty DataFrame to store the cumulative 5th columns
    file_path = folder_path

    df = pd.read_csv(file_path, sep=" ")  # Adjust the separator if needed

    observatory_column = df.columns[-3]
    df[WhatToPlot] = df[WhatToPlot] * 100
    grouped_data = df.groupby(observatory_column)

    time_shifted = pd.to_datetime(df["Time[UTC]"]) + pd.to_timedelta(
        df["Duration"], unit="m"
    )
    df["end_obs"] = time_shifted
    sorted_df = df.sort_values(by="end_obs")

    sorted_df["cumsum"] = sorted_df[
        WhatToPlot
    ].cumsum()  # Cumulative sum of PGal values
    max_value_indices = sorted_df.groupby("end_obs")["cumsum"].idxmax()
    max_rows = sorted_df.loc[max_value_indices]
    # Dictionary to track if legend has been added for each observatory
    legend_added = {}
    for index, (observatory, observatory_data) in enumerate(grouped_data):
        duration = observatory_data["Duration"].iloc[0]

        # Extract relevant columns for the plot
        x_values = pd.to_datetime(observatory_data["Time[UTC]"])
        y_values = observatory_data[WhatToPlot].cumsum()

        #  Calculate shifted x_values
        y_values_shifted = [0] + y_values.tolist()[:-1]
        x_shifted = x_values + pd.Timedelta(minutes=duration)

        # Initialize the legend flag for the observatory if not already initialized
        if observatory not in legend_added:
            legend_added[observatory] = False

        # Draw lines between corresponding points
        for xv, yvs, xs, yv in zip(x_values, y_values_shifted, x_shifted, y_values):
            if not legend_added[observatory]:
                plt.plot(
                    [xv, xs],
                    [yvs, yv],
                    marker=".",
                    linestyle="-",
                    color=colors[index],
                    label=observatory,
                )
                legend_added[observatory] = True
            else:
                plt.plot(
                    [xv, xs], [yvs, yv], marker=".", linestyle="-", color=colors[index]
                )

    plt.plot(
        max_rows["end_obs"],
        max_rows["cumsum"],
        marker="+",
        linestyle="-",
        label="All",
        color="black",
    )
    plt.gca().xaxis.set_major_locator(mdates.AutoDateLocator())
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
    plt.gcf().autofmt_xdate()  # Auto-format the date labels for better readability
    plt.xlabel("Date [UTC]")
    plt.ylabel("Total Covered Probability [%]")
    plt.xticks(rotation=45)  # Rotate x-axis labels for better readability
    plt.legend()

    plt.grid(which="both")
    if interactive:
        plt.show()
    else:
        plt.savefig(
            "%s_%s_cumulativeDistribution.png" % (event_name, WhatToPlot),
            dpi=300,
            bbox_inches="tight",
        )
        plt.close()


def ObtainCumulativeProbabilityPlotLog(
    folder_path, event_name, WhatToPlot, interactive=False
):

    # Initialize an empty DataFrame to store the cumulative 5th columns
    file_path = folder_path

    df = pd.read_csv(file_path, sep=" ")  # Adjust the separator if needed

    observatory_column = df.columns[-3]
    df[WhatToPlot] = df[WhatToPlot] * 100
    grouped_data = df.groupby(observatory_column)

    time_shifted = pd.to_datetime(df["Time[UTC]"]) + pd.to_timedelta(
        df["Duration"], unit="m"
    )
    df["end_obs"] = time_shifted
    sorted_df = df.sort_values(by="end_obs")

    sorted_df["cumsum"] = sorted_df[
        WhatToPlot
    ].cumsum()  # Cumulative sum of PGal values
    max_value_indices = sorted_df.groupby("end_obs")["cumsum"].idxmax()
    max_rows = sorted_df.loc[max_value_indices]
    # Dictionary to track if legend has been added for each observatory
    legend_added = {}
    for index, (observatory, observatory_data) in enumerate(grouped_data):
        duration = observatory_data["Duration"].iloc[0]

        # Extract relevant columns for the plot
        x_values = pd.to_datetime(observatory_data["Time[UTC]"])
        y_values = observatory_data[WhatToPlot].cumsum()

        #  Calculate shifted x_values
        y_values_shifted = [0] + y_values.tolist()[:-1]
        x_shifted = x_values + pd.Timedelta(minutes=duration)

        # Initialize the legend flag for the observatory if not already initialized
        if observatory not in legend_added:
            legend_added[observatory] = False

        # Draw lines between corresponding points, skipping the first point
        for i, (xv, yvs, xs, yv) in enumerate(
            zip(x_values, y_values_shifted, x_shifted, y_values)
        ):
            if i == 0:
                continue  # Skip the first point
            if not legend_added[observatory]:
                plt.plot(
                    [xv, xs],
                    [yvs, yv],
                    marker=".",
                    linestyle="-",
                    color=colors[index],
                    label=observatory,
                )
                legend_added[observatory] = True
            else:
                plt.plot(
                    [xv, xs], [yvs, yv], marker=".", linestyle="-", color=colors[index]
                )

    plt.plot(
        max_rows["end_obs"],
        max_rows["cumsum"],
        marker="+",
        linestyle="-",
        label="All",
        color="black",
    )
    plt.gca().xaxis.set_major_locator(mdates.AutoDateLocator())
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
    plt.gcf().autofmt_xdate()  # Auto-format the date labels for better readability
    plt.xlabel("Date [UTC]")
    plt.ylabel("Total Covered Probability [%]")
    plt.xticks(rotation=45)  # Rotate x-axis labels for better readability
    plt.legend()
    plt.yscale("log")  # Set Y-axis to log scale
    plt.grid(which="both")
    if interactive:
        plt.show()
    else:
        plt.savefig(
            "%s_%s_cumulativeDistribution.png" % (event_name, WhatToPlot),
            dpi=300,
            bbox_inches="tight",
        )
        plt.close()
