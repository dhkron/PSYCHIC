function [a_b,a_t,a_c] = GetTADsMask(raw_domains,MIN_DIAG,MAX_DIAG,a_size,res)

domains = floor(raw_domains/res)+1; %To start from 1
domains_map_xor = zeros(a_size);
domains_map_one = zeros(a_size);
for i = MIN_DIAG:MAX_DIAG
	for i = MIN_DIAG:MAX_DIAG
		start_pos = a_size*i + 1;
		domains_diag_xor = zeros(1,a_size-i);
		domains_diag_one = zeros(1,a_size-i);
		for domain = domains'
			domain(2) = min(domain(2),a_size); %Fix overflow

			bnd_start = domain(1);
			bnd_end = domain(2)-i;
			min_bnd = min(bnd_start,bnd_end);
			min_bnd = max(1, min_bnd);
			max_bnd = max(bnd_start,bnd_end);
			max_bnd = min(max_bnd, a_size-i);
			domains_diag_xor(min_bnd:max_bnd) = ~domains_diag_xor(min_bnd:max_bnd);

			min_bnd = max(1, bnd_start);
			max_bnd = min(bnd_end, a_size-i);
			domains_diag_one(min_bnd:max_bnd) = 1; %Do not extend TADs to tiles

			%For calculating TAD strength vs. Length
		end
		domains_map_xor(start_pos:a_size+1:end) = domains_diag_xor;
		domains_map_one(start_pos:a_size+1:end) = domains_diag_one;
		%domains_map(start_pos:a_size+1:end) = round(rand(1,a_size-i)); %Random model
	end
	a_c = domains_map_xor-domains_map_one; %Just TAD interactions, not TADs
	a_t = domains_map_one;
	a_b = 1-domains_map_xor;
end

end
