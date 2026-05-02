from PyQt6.QtWidgets import QMainWindow
from sequence_analysis_gui import Ui_Sequence_analysis_window
from RNA_protein_Class import RNASequence


class LogicController:
    """
    Controls interaction between the input window and analysis window.

    Responsibilities:
    - Handle user input from the RNA input window
    - Validate and create RNASequence objects
    - Open and manage the analysis window
    - Execute selected analysis methods
    - Display results in the GUI
    """

    def __init__(self, ui, window) -> None:
        """
        Initialize controller.

        Args:
            ui: Input window UI object
            window (QMainWindow): Input window instance
        """
        self.ui = ui
        self.window: QMainWindow = window
        self.rna_sequence: RNASequence | None = None

        self.analysis_window: QMainWindow | None = None
        self.analysis_ui: Ui_Sequence_analysis_window | None = None

        self.connect_signals()

    def connect_signals(self) -> None:
        """
        Connect input window buttons to controller methods.
        """
        self.ui.analyze_button.clicked.connect(self.open_analysis_window)
        self.ui.pushButton.clicked.connect(self.clear_input)

    def clear_input(self) -> None:
        """
        Clear the RNA input text field.
        """
        self.ui.RNA_input.clear()

    def open_analysis_window(self) -> None:
        """
        Validate RNA input and open the analysis window.

        Error displays if invalid input is given in the input box.
        """
        sequence: str = self.ui.RNA_input.toPlainText().strip()

        try:
            self.rna_sequence = RNASequence(sequence)
        except ValueError:
            self.ui.RNA_input.setText("Invalid RNA sequence.")
            return

        self.analysis_window = QMainWindow()
        self.analysis_ui = Ui_Sequence_analysis_window()
        self.analysis_ui.setupUi(self.analysis_window)

        self.analysis_ui.METHOD_RESULTS_txt.setText("")
        self.analysis_ui.run_methods_button.clicked.connect(self.run_selected_methods)
        self.analysis_ui.return_button.clicked.connect(self.return_to_input)

        self.analysis_window.show()
        self.window.hide()

    def return_to_input(self) -> None:
        """
        Close analysis window and return to input window.
        """
        if self.analysis_window:
            self.analysis_window.close()
        self.window.show()

    def run_selected_methods(self) -> None:
        """
        Execute selected analysis methods and display results.

        Builds a formatted results string based on selected checkboxes.
        """
        if not self.rna_sequence or not self.analysis_ui:
            return

        results: str = ""

        if self.analysis_ui.show_codon_checkbox.isChecked():
            codons: list[str] = self.rna_sequence.get_codons()
            results += "Codons:\n"
            results += " | ".join(codons) + "\n\n"

        if self.analysis_ui.translate_codon_checkbox.isChecked():
            protein: list[str] = self.rna_sequence.translate()
            results += "Translation:\n"
            results += "-".join(protein) + "\n\n"

        if self.analysis_ui.find_start_codon_checkbox.isChecked():
            start_codon: int = self.rna_sequence.get_start_codon()
            results += "Start Codon:\n"

            if start_codon == -1:
                results += "No Start Codon found.\n\n"
            else:
                results += f"AUG found at index {start_codon}.\n\n"

        if self.analysis_ui.amino_acid_count_checkbox.isChecked():
            amino_count: dict[str, int] = self.rna_sequence.amino_acids_count()
            results += "Amino Acid Counts:\n"

            if not amino_count:
                results += "No amino acids detected.\n\n"
            else:
                for amino, count in amino_count.items():
                    results += f"{amino}: {count}\n"
                results += "\n"

        if self.analysis_ui.product_prediction_checkbox.isChecked():
            prediction: str = self.rna_sequence.predict_protein_functionality()
            results += "Product Prediction:\n"
            results += f"{prediction}\n\n"

        if results == "":
            results = "Please select at least one analysis method."

        self.analysis_ui.METHOD_RESULTS_txt.setText(results)
