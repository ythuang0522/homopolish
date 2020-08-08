import os

class FileManager:
    """
    Performs simple OS operations like creating directory for output or path to all the files.
    """
    @staticmethod
    def handle_output_directory(output_dir):
        """
        Process the output directory and return a valid directory where we save the output
        :param output_dir: Output directory path
        :return:
        """
        # process the output directory
        if output_dir[-1] != "/":
            output_dir += "/"
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        return os.path.abspath(output_dir)
