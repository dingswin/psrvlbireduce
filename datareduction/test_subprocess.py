#!/usr/bin/env python
import subprocess
difmap = subprocess.Popen('difmap', stdin=subprocess.PIPE, text=True)
difmap.stdin.write("float pkflux\n")
difmap.stdin.write("float peakx\n")
difmap.stdin.write("exit\n")
#difmap.communicate(input="float pkflux\nfloat peakx\nexit\n")
#difmap.poll()
difmap.stdin.close()
#if difmap.wait() != 0:
#    print("There were some errors")
difmap.wait()
#difmap.terminate()
print('checking')
