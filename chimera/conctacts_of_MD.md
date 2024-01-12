Steps:
  1. Create selections in Chimera GUI by selecting the first group and then introducing the following ``namesel Group1`` command. Do the same with the second group but this time use the ``namesel Group2`` command
  2. Go to the `MD Movie` window, then to the `Per-frame` panel and choose `Define script`. Change the `Interpret script as` to Python.
  3. Paste the following script and introduce a route and name (like route_to_dir/filename) of where the outputs should be written. It will be run in each step.
  ``python
from chimera import runCommand

route_file = "" # Change the route!! Input the route from the root.

runCommand("findclash FLAP test ALOX5 overlapCutoff -0.4 hbondAllowance 0 bondSeparation 4 intraMol false saveFile " + route_file + str(mdInfo['frame']) + '.dat')
``

  4. Play the MD simulation. You can choose to run a part of the MD and/or to use steps to skip some frames. The script will be only executed for each frame that is represented.
