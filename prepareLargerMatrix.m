function [input_data, sim_data] = prepareLargerMatrix(final_cluster_genes_data, samples, num_clusters, num_samples, simulated_data)
    input_data = [];
    sim_data = [];

    for i4 = 1:num_samples
        eval(['A = final_cluster_genes_data{samples(i4)};'])
        A = A(~cellfun(@isempty,A));
        tmp = [];

        for j4 = 1:num_clusters
            aa = A{j4};
            gene_data = cell2mat(aa(:, 2:end));
            time_point_means = mean(gene_data, 1);
            correlations = zeros(size(aa, 1), 1);

            for k = 1:size(aa, 1)
                correlations(k) = corr(gene_data(k, :)', time_point_means');
            end

            [~, idx1] = sort(correlations, 'descend');

            if length(idx1) <= 5
                random_pick1 = idx1(randi(length(idx1)));
                tmp = [tmp; gene_data(random_pick1, :)];
            else
                random_pick1 = idx1(randi(5));
                tmp = [tmp; gene_data(random_pick1, :)];
            end
        end

        input_data = [input_data; tmp];
        eval(['B = simulated_data{samples(i4)};'])
        sim_data = [sim_data; B];
    end
end
