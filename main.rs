use rand::Rng;
use scoped_threadpool;
use scoped_threadpool::Pool;
use num_cpus;

fn main()
{
    let mut pool = Pool::new(num_cpus::get() as u32);
    let dimensions = 10;
    
    pool.scoped(|scope| {
        scope.execute(move || {
            let mut pluszle_not_found = true;
            let mut number_twos = 0;
            let mut dif: f32 = 0.0;
            let mut loop_counter = 0;
            let mut all_solutions;
            let mut solution:[Vec<Vec<Vec<u32>>>; 2] = [Vec::new(), Vec::new()];
            let mut final_pluszle   = (Vec::new(), Vec::new());
        while pluszle_not_found
        {
            dif = 0.0;
            number_twos = 0;
            let (pluszle_matrix, pluszle_sums) = make_sums(dimensions, dimensions);
            let pluszle_rows = row_vectors(&pluszle_matrix, dimensions);
            let pluszle_cols = col_vectors(&pluszle_matrix, dimensions);
            let pluszle = [pluszle_rows, pluszle_cols];

            //println!("{:?} {:?} {:?}", pluszle_matrix, pluszle_sums, pluszle[0]);
            all_solutions = [Vec::new(), Vec::new()];

            for n in 0..all_solutions.len()
            {
                for (index, rc) in pluszle[n].iter().enumerate()
                {
                    let sol = get_solutions(rc, pluszle_sums[index+(n*dimensions)]);
                    all_solutions[n].push(sol);
                }
            }

            // Create global masks
            let mut found_solution = false;
            loop_counter = 0;

            while !found_solution && loop_counter < 20
            {
                let mut masks = [Vec::new(), Vec::new()];
                for n in 0..2
                {
                    for sol in &all_solutions[n]
                    {
                        let (mask_temp, num) = create_global_mask(sol, dimensions);
                        if loop_counter == 0
                        {
                            number_twos += num;
                        }
                        masks[n].push(mask_temp);
                    }
                }

                // Change masks
                let mut changed = false;
                for x in 0..dimensions
                {
                    for y in 0..dimensions
                    {
                        let a = redefine_orthognal_masks(&mut masks, y, x);
                        match a {
                            Some((h, v)) => {
                                if masks[0][x] != h || masks[1][y] != v
                                {
                                    changed = true;
                                    dif += 1.0
                                }
                                masks[0][x] = h;
                                masks[1][y] = v;
                            },
                            None => println!("Error")
                        }
                    }
                }

                if !changed
                {
                    loop_counter = 11000;
                }

                if found_solution == false
                {
                    // Remove solutions which do not fit the new mask
                    let mut newest_solutions = [Vec::new(), Vec::new()];
                    found_solution = true;
                    for n in 0..2
                    {
                        for (index, sol) in all_solutions[n].iter().cloned().enumerate()
                        {
                            let mut new_solutions = Vec::new();
                            for rc in sol
                            {
                                if masks[n][index].iter().zip(&rc).filter(|&(a, b)| a != b && *a != 2_u32).count() > 0
                                {
                                    continue;
                                } else {
                                    let mask_sum = rc.iter().cloned().sum::<u32>();
                                    let row_len = rc.len() as u32;
                                    if (mask_sum != 0 && mask_sum != 1 && mask_sum != row_len && mask_sum != (row_len - 1)) || dimensions < 4
                                    {
                                        new_solutions.push(rc);
                                    } else {
                                        new_solutions.push(rc);
                                        found_solution = false;
                                    }
                                }
                                if new_solutions.len() > 1
                                {
                                    found_solution = false;
                                }
                            }

                            newest_solutions[n].push(new_solutions);
                        }
                        all_solutions[n] = newest_solutions[n].clone();
                    }
                }
                loop_counter += 1;
            }
            
            if loop_counter < 20
            {
                pluszle_not_found = false;
                solution = all_solutions;
                final_pluszle = (pluszle_matrix, pluszle_sums);
            }

        }
            println!("{:?} {:?} {:?} {:?} {:?} {:?}", (dif * ((number_twos as f32) / 2 as f32)) / (dimensions.pow(4) as f32), dif, number_twos, loop_counter, final_pluszle, solution);

        });
    });
}

fn row_vectors(matrix:&Vec<u32>, num_rows: usize) -> Vec<Vec<u32>>
{
    let mut return_vector: Vec<Vec<u32>> = Vec::new();
    for n in (0..matrix.len()).step_by(num_rows)
    {
        let r: Vec<u32>  = matrix[n..n+num_rows].to_vec();
        return_vector.push(r);
    }

    return_vector
}

fn col_vectors(matrix:&Vec<u32>, num_cols: usize) -> Vec<Vec<u32>>
{
    let mut return_vector: Vec<Vec<u32>> = Vec::new();
    for n in 0..num_cols
    {
        let c: Vec<_> = matrix.iter().cloned().skip(n).step_by(num_cols).collect();
        return_vector.push(c);
    }
    return return_vector;
}

fn redefine_orthognal_masks(masks: &mut[Vec<Vec<u32>>;2], x: usize, y: usize) -> Option<(Vec<u32>, Vec<u32>)>
{
    let hor_mask= &mut masks[0][y].clone();
    let ver_mask = &mut masks[1][x].clone();
    if hor_mask[x] != ver_mask[y] && (hor_mask[x] == 2 || ver_mask[y] == 2)
    {
        if hor_mask[x] != 2
        {
            ver_mask[y] = hor_mask[x];
        }
        else if ver_mask[y] != 2
        {
            hor_mask[x] = ver_mask[y];
        }

        Some((hor_mask.to_vec(), ver_mask.to_vec()))

    }
    else if hor_mask[x] == ver_mask[y]
    {
        Some((hor_mask.to_vec(),ver_mask.to_vec()))
    }
    else {
        None
    }


}

fn create_global_mask(solution_mask: &Vec<Vec<u32>>, mask_len: usize)   -> (Vec<u32>, i32)
{
    let mut return_mask = Vec::new();
    let mut new_mask  = Vec::new();
    let mut number_twos = 0;
    for n in 0..mask_len
    {
        let mut inverse_vec = Vec::new();
        for mask in solution_mask
        {
            inverse_vec.push(mask[n]);
        }
        new_mask.push(inverse_vec);
    }
    for nm in new_mask
    {
        if nm.windows(2).all(|x| x[0] == x[1])
        {
             return_mask.push(nm[0]);
        }
        else
        {
            return_mask.push(2);
            number_twos += 1;
        }
    }
    (return_mask, number_twos)
}


fn get_solutions(line_vec: &Vec<u32>, sum: u32) -> Vec<Vec<u32>>
{
    let len_ = line_vec.len();
    let masks: Vec<Vec<u32>> = bin_mask(len_);
    let mut solutions: Vec<Vec<u32>> =  Vec::new();
    for mask in masks
    {
        let mut part_sum : Vec<u32> = Vec::new();
        for m in 0..mask.len()
        {
            if mask[m] == 1
            {
               part_sum.push(line_vec[m])
            }
        }
        if part_sum.iter().sum::<u32>() == sum
        {
            solutions.push(mask);
        }
    }
    solutions
}

fn bin_mask(size: usize)  -> Vec<Vec<u32>>
{
    let mut masks = Vec::new();
    for number in 1..2_i32.pow(size as u32)
    {
        let bin_string = format!("{:0in_size$b}", number, in_size = size);
        let bin_vec: Vec<u32> = bin_string.chars().map(|c| c.to_digit(2).unwrap()).collect();
        masks.push(bin_vec);

    }
    masks
}

fn make_sums(rows: usize, cols: usize) -> (Vec<u32>, Vec<u32>)
{
    let mut false_vector = true;
    let mut sums: Option<Vec<u32>>;
    let mut sum_new = Vec::new();
    let mut new_matrix = Vec::new();

    while false_vector
    {
        let mut matrix: Vec<u32> = Vec::new();
        let mut pattern: Vec<u32> = Vec::new();
        let mut rng = rand::thread_rng();

        for _ in 0..cols * rows
        {
            matrix.push(rng.gen_range(1, 10));
            pattern.push(rng.gen_range(0, 2));
        }


        sums = get_sums_vectors(&matrix, &pattern, cols);
      
        match sums {
            Some(s) => {false_vector = false; sum_new = s; new_matrix = matrix},
            None => {false_vector = true}
        }
    }

    (new_matrix, sum_new)

}


fn get_sums_vectors(matrix: &Vec<u32>, pattern: &Vec<u32>, num_cols: usize) -> Option<Vec<u32>>
{
    let iter = pattern.iter().zip(matrix.iter());
    let new_matrix: Vec<u32> = iter.map(|(&x, &y)| x*y).collect();
    let mut row_sums: Vec<u32> = Vec::new();
    let mut col_sums: Vec<u32> = Vec::new();

    for sl in new_matrix.chunks(num_cols)
    {

        row_sums.push(sl.iter().sum());
    }

    for col in 0..num_cols
    {
        col_sums.push(new_matrix.iter().skip(col).step_by(num_cols).sum());
    }

     row_sums.append(&mut col_sums);


     if row_sums.contains(&(0 as u32))
     {
         None
     }
    else {
        Some(row_sums)
    }



}