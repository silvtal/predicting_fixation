def generate_PCGtable(N, Av, L=None, LN=None, outfile=None, seed=None):
    """
    Generate a table with N rows, where each row represents a group with an
    average value and a number of leaves.

    Args:
        N (int): Number of groups.
        Av (list of floats): List of niche sizes (size for each group). The length
            of the list must be equal to N, and the values must sum to 1.
        L (list of ints, optional): List of number of leaves for each group. If this
            argument is provided, it takes precedence over the LN argument. Default value
            is None.
        LN (int, optional): Total number of leaves. If this argument is provided,
            the number of leaves for each group will be randomly generated. Default
            value is None.
        outfile (str, optional): Output file path. If provided, the generated table
            will be saved in CSV format at the specified location. Default value is None.
        seed (int, optional): Random seed. Default value is None.

    Returns:
        list of lists: A list of N+1 lists, where each list represents a row in the
        generated table. The first two elements of each row are the group name and the
        average value, respectively. The third element is either the number of leaves
        for the group or a semicolon-separated string of leaf names. The last row of the
        table contains the total number of leaves.

    Raises:
        ValueError: If both L and LN arguments are None.
    
    If LN argument is None, the function randomly generates the number of leaves for
    each group. The total number of leaves is determined by the L argument. The function
    ensures that each group has at least one leaf and the total number of leaves is
    distributed evenly among the groups. The number of leaves for the last group is
    the remainder of the division, and is added to the last row of the table as the
    total number of leaves.
    
    Example usage:
	N = 3
	Av = [0.2, 0.3, 0.5]
	LN = 10

	table = generate_table(N, Av, LN=LN)

	for row in table:
	    print("\t".join(str(x) for x in row))
    
    Example usage without LN argument:
	N = 4
	Av = [0.1, 0.2, 0.3, 0.4]
	L = [1, 2, 3, 4]

	table = generate_table(N, Av, L=L)

	for row in table:
	    print("\t".join(str(x) for x in row))

    """
    
    import random
    import csv
    
    if seed is not None:
        random.seed(seed)
    if L is None:
        if LN is None:
            raise ValueError("Either L or LN argument must be provided.")
        else:
            L = [LN // N] * N
            if sum(L) < LN:
                L[-1] += LN - sum(L) # do not leave out any leaves
    else:
        if LN is not None and sum(L) != LN:
            print("Warning: L argument takes precedence over LN argument.")
        LN = sum(L)

    # Create list of leaves
    leaves = [f"otu{i}" for i in range(1, LN+1)]

    # Create table
    table = []
    for i in range(N):
        # Select leaves for this row
        row_leaves = random.sample(leaves, L[i]) # if we randomly pick from a 
        for leaf in row_leaves:                  # lognorm, we expect the groups
            leaves.remove(leaf)                  # to have a lognorm distrb too
        # Append row to table
        row = ["group" + str(i+1), Av[i], ";".join(row_leaves)]
        table.append(row)
    
    table.append([LN, "", ""])

    if outfile is not None:
        with open(outfile, "w", newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['Core', 'Average', 'Leaves'])
            writer.writerows(table)
    return table
