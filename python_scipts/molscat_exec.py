from subprocess import PIPE, Popen
from time import sleep
import time
from python_scipts.input_creator import Input


class CygwinExecutor:
    path = "./projects/molscat/"

    def execute(self, input: Input, filename: str) -> str:
        p = Popen("C:/cygwin64/bin/bash.exe --login", stdin=PIPE, stdout=PIPE)
        
        input_filepath = f"inputs/{filename}.input"
        output_filepath = f"outputs/{filename}.output"

        with open(input_filepath, "w") as file:
            file.write(input.to_string())

        p.stdin.write(f"cd {self.path}\n".encode()) # type: ignore
        p.stdin.write(f"time source_code/molscat-basic < {input_filepath} > {output_filepath}".encode()) # type: ignore
        p.stdin.close() # type: ignore
        p.wait()

        out = p.stdout.read() # type: ignore
        p.terminate()

        return out.decode()