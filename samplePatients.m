function samples = samplePatients(random_patients, num_samples)
    elements = random_patients;
    num_repeats = floor(num_samples / numel(elements));
    samples = zeros(1, num_samples);
    samples(1:numel(elements) * num_repeats) = repmat(elements, 1, num_repeats);
    remaining_samples = num_samples - numel(elements) * num_repeats;

    if remaining_samples > 0
        samples(end - remaining_samples + 1:end) = elements(randperm(numel(elements), remaining_samples));
    end

    samples = samples(randperm(num_samples));
end
