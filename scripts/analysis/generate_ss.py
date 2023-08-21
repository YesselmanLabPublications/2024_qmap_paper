import subprocess
import pandas as pd
import click


class Varna(object):
    def __init__(self):
        self.varna_command = (
            "java -cp ./resources/VARNA.jar fr.orsay.lri.varna.applications.VARNAcmd"
        )

    def new_image(
        self,
        filename,
        sequence,
        structure,
        colormapstyle="red",
        highlight_res=[],
        highight_bp=[],
    ):
        colormap = [0 for x in range(len(sequence))]
        colormap_str = ";".join([str(x) for x in colormap])

        subprocess.call(
            self.varna_command
            + f" -sequenceDBN '{sequence}'"
            + f" -structureDBN '{structure}'"
            + f" -o {filename}",
            shell=True,
        )


@click.command()
@click.argument("csv_file", type=click.Path(exists=True))
def main(csv_file):
    """
    main function for script
    """
    df = pd.read_csv(csv_file)
    v = Varna()
    for i, row in df.iterrows():
        v.new_image(
            f"resources/images/{row['name']}.png", row["act_seq"], row["act_ss"]
        )


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    main()
