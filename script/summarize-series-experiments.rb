# Script to generate table
gem 'parallel'
gem 'descriptive_statistics'
require 'parallel'
require 'descriptive_statistics'

results_dir = File.expand_path(ARGV[0])
n_rounds = 10
n_forks = (10..100).step(5)
periods = [
  0,
  100,
  1000,
  10_000,
  100_000
]

p 'number of forks, steps per iteration, total iterations, average inter origin distance'
Parallel.each(n_forks) do |f|
  Parallel.each(periods) do |period|
    p "Forks: #{f} | Period: #{period}"
    v_iterations = []
    v_iod = []
    (0..n_rounds - 1).each { |round|
      (0..1000 - 1).each do |sim|
        begin
          cell_file = "#{results_dir}/round_#{round}_false_#{f}_#{period}/simulation_#{sim}/cell.txt"
          (iterations, iod) =  File.open(cell_file).each("\t").map(&:to_i).to_a.slice(2, 2)
          v_iterations.append(iterations)
          v_iod.append(iod)
        rescue Errno::ENOENT => e
          # p "#{cell_file} not found"
        rescue e
          p e
        end
      end
    }
    p "#{f}, #{period}, #{v_iterations.mean}, #{v_iterations.standard_deviation}, #{v_iod.mean}, #{v_iod.standard_deviation}"
  end
end
