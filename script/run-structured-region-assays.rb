#!/usr/bin/ruby

require 'optparse'
require 'ruby-progressbar'

options = {
  simulations: 1000,
  timeout: 100_000_000,
  speed: 1,
  threads: 20,
  organism: 'TcruziCLBrenerEsmeraldo-like'
}

OptionParser.new do |opts|
  opts.banner = 'Usage: simulations.rb [options]'

  opts.on('-n', '--cells [Integer]', 'Number of simulations to run') do |v|
    options[:simulations] = Integer(v)
  end

  opts.on('-t', '--timeout [Integer]', 'Simulation timeout') do |v|
    options[:timeout] = Integer(v)
  end

  opts.on('-p', '--threads [Integer]', 'Number of threads to use') do |v|
    options[:threads] = Integer(v)
  end

  opts.on('-o', '--organism [Integer]', 'Organism name') do |v|
    options[:organism] = v
  end

  opts.on('-s', '--speed [Integer]', 'Replisome speed') do |v|
    options[:speed] = Integer(v)
  end

  opts.on('-r', '--round [Integer]', 'The number of the round in a series of similar experiments') do |v|
    options[:round] = Integer(v)
  end

  opts.on('-h', '--help', 'Prints this help') do
    puts opts
    exit
  end
end.parse!

replisome_count = (10..100).step(10)

transcription_period = [
  0,
  100,
  1000,
  10_000,
  100_000
]

pb = ProgressBar.create(
  title: "Completed Simulations", 
  total: replisome_count.size * transcription_period.size,
  format: '%t %c/%C |%B| %p%% %a %E'
)


pb.log "Starting simulations for round #{options[:round]}\n"
sleep(3)
start_time = Process.clock_gettime(Process::CLOCK_MONOTONIC)
transcription_period.each do |period|
  replisome_count.each do |replisomes|
    pb.log "Running for replisomes=#{replisomes} and period=#{period}\n"
    Dir.chdir('./build') do
      system("mkdir -p output_structured_regions/round_#{options[:round]}_false_#{replisomes}_#{period}")
      system("nice -n 20 ./simulator --cells #{options[:simulations]} \
        --organism '#{options[:organism]}' \
        --resources #{replisomes} \
        --speed #{options[:speed]} \
        --period #{period} \
        --timeout #{options[:timeout]} \
        --threads #{options[:threads]} \
        --name round_#{options[:round]} \
        --output output_structured_regions \
        1> output_structured_regions/round_#{options[:round]}_false_#{replisomes}_#{period}/simulation_out \
        2> output_structured_regions/round_#{options[:round]}_false_#{replisomes}_#{period}/simulation_err \
        ")
      pb.increment
    end
  end
end

end_time = Process.clock_gettime(Process::CLOCK_MONOTONIC)
elapsed_time = end_time - start_time
File.open('./build/output_structured_regions/elapsed_time', 'w') do |file|
  file.write("Elapsed time in seconds: #{elapsed_time}\n")
end
