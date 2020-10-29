#!/usr/bin/python3

import argparse, os, subprocess, shutil


def check_runtime(l, t, d, path, time_bd, repeat):
  bench = "bench.csv"
  warmup = 10
  time_bd_total = 2 * time_bd * (warmup + repeat)

  newton_cmd = "./newton -l{l} -t{t} {d}".format(l = l, t = t, d = d)
  cmd = "hyperfine --export-csv {bench} --time-unit millisecond --warmup {warmup} --max-runs {repeat} \"{newton_cmd}\"".format(bench = "bench.csv", warmup = warmup, repeat = repeat, newton_cmd = newton_cmd)
  subprocess.Popen(cmd, shell = True, cwd=path, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL).wait(time_bd_total)
 
  try:
    with open(os.path.join(path, bench), "r") as f:
      time = float(f.readlines()[-1].split(",")[-2])
      # print("CALIBRATE: {} -> {}".format(newton_cmd, 2*time))
      if time > time_bd:
        print("TOO SLOW ({} ms) FOR {}".format(time, newton_cmd))
        return False
      else:
        return True
  except:
    pass

  print("RUNTIME ERROR FOR {}".format(newton_cmd))
  return False

def check_runtimes(path):
  return (   check_runtime(1000, 1, 1, path, 0.04, 100)
         and check_runtime(1000, 1, 2, path, 0.11, 100)
         and check_runtime(1000, 1, 5, path, 0.44, 100)
         and check_runtime(1000, 1, 7, path, 0.80, 100)

         and check_runtime(1000, 1, 5, path, 0.44, 100)
         and check_runtime(1000, 2, 5, path, 0.23, 100)
         and check_runtime(1000, 3, 5, path, 0.16, 100)
         and check_runtime(1000, 4, 5, path, 0.12, 100)

         and check_runtime(1000, 10, 7, path, 0.10, 100)
         and check_runtime(2000, 10, 7, path, 0.36, 100)
         and check_runtime(3000, 10, 7, path, 0.80, 100)
         )

def run_check_build():
  # check files
  is_valid_file = lambda f: ( 
       f in ["makefile", "Makefile"] or
       f.endswith(".cc") or f.endswith(".c") or
       f.endswith(".hh") or f.endswith(".h") )
  is_valid_folder = lambda folder: (
    all( all(map(is_valid_file, files))
         for (root, _, files) in os.walk(folder) ) )

  print( "checking for additional files..." )
  if not is_valid_folder(extraction_path):
    print("ADDITIONAL FILES IN TAR")
    return False

  # build clean build
  print( "bulding and cleaning..." )
  subprocess.Popen(["make", "newton"], cwd=extraction_path).wait()
  subprocess.Popen(["make", "clean"], cwd=extraction_path).wait()
  if not is_valid_folder(extraction_path):
    print("ADDITIONAL FILES AFTER BUILD CLEAN")
    return False

  print( "bulding..." )
  subprocess.Popen(["make", "newton"], cwd=extraction_path, stdout=subprocess.DEVNULL).wait()
  return True

def run_copy_picture():
  print( "creating test picture..." )
  subprocess.run([extraction_path + "/newton", "-l1000", "-t4", "7"])
  subprocess.run(["convert", "newton_attractors_x7.ppm", "newton_attractors_x7.png"])
  subprocess.run(["convert", "newton_convergence_x7.ppm", "newton_convergence_x7.png"])
  subprocess.run(["mv", "newton_attractors_x7.png", "pictures/" + stem + "_attractors.png"])
  subprocess.run(["mv", "newton_convergence_x7.png", "pictures/" + stem + "_convergence.png"])
  subprocess.run(["rm", "newton_attractors_x7.ppm"])
  subprocess.run(["rm", "newton_convergence_x7.ppm"])

def run_check_runtime():
  print( "checking runtimes..." )
  return check_runtimes(extraction_path)



subprocess.run(["mkdir", "-p", "extracted"])
subprocess.run(["mkdir", "-p", "pictures"])

parser = argparse.ArgumentParser()
parser.add_argument("tarfile")
args = parser.parse_args()

tar_file = args.tarfile
assert( tar_file.endswith(".tar.gz") )

stem = tar_file[:-7]
while stem.find("/") != -1:
  stem = stem[stem.find("/")+1:]

extraction_path = "extracted/" + stem
if os.path.isdir(extraction_path):
  print("FATAL: extreaction path {} exists already".format(extraction_path))
  exit(1)

print( "extracting..." )
subprocess.run(["mkdir", extraction_path])
subprocess.run(["tar", "xf", tar_file, "-C", extraction_path])

passed = run_check_build()
print(f"Passed check build: {passed}")
passed = run_copy_picture()
print(f"Passed copy picture: {passed}")
passed = run_check_runtime()
print(f"Passed run times: {passed}")

# clean
print( "final cleaning..." )
shutil.rmtree(extraction_path)

# feedback summary
if passed:
  print("submission passes script")
else:
  print("submission DOES NOT pass script")
