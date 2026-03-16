from __future__ import annotations

import unittest
from unittest.mock import Mock, patch

from screened_bond_valence.cli import main


class CliTests(unittest.TestCase):
    def test_batch_command_runs_processor_for_all_pairs(self) -> None:
        processor = Mock()

        with patch("screened_bond_valence.cli._build_processor", return_value=processor) as processor_factory:
            exit_code = main(
                [
                    "batch",
                    "--api-key",
                    "test-key",
                    "--cations",
                    "Li",
                    "Na",
                    "--anions",
                    "O",
                    "--algorithms",
                    "shgo",
                ]
            )

        self.assertEqual(exit_code, 0)
        processor_factory.assert_called_once()
        self.assertEqual(processor.process_cation_system.call_count, 2)
        processor.process_cation_system.assert_any_call("Li", "O")
        processor.process_cation_system.assert_any_call("Na", "O")

    def test_web_command_launches_app(self) -> None:
        with patch("screened_bond_valence.cli._launch_web_app") as launch_app:
            exit_code = main(
                [
                    "web",
                    "--server-name",
                    "0.0.0.0",
                    "--server-port",
                    "9000",
                    "--share",
                ]
            )

        self.assertEqual(exit_code, 0)
        launch_app.assert_called_once_with(
            server_name="0.0.0.0",
            server_port=9000,
            share=True,
            inbrowser=False,
        )


if __name__ == "__main__":
    unittest.main()
