from subprocess import Popen, PIPE
import os

path = "./projects/molscat/"
p = Popen("C:/cygwin64/bin/bash.exe --login", stdin=PIPE, stdout=PIPE)

p.stdin.write(f"cd {path}\n".encode()) # type: ignore
p.stdin.write(b"source_code/molscat-basic < examples/input/molscat-2chan-LJ.input > outputs/2chan.output") # type: ignore
p.stdin.close() # type: ignore
p.wait()
out = p.stdout.read() # type: ignore
print(out.decode())
p.terminate()