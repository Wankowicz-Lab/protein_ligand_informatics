

const matrix_width = 4;
const job_array = matrix_width * (matrix_width - 1) / 2 // 1-6

// matrix is mirrored along the diagonal
// do not compare same index with itself
// 1-2, 1-3, 1-4, 2-3, 2-4, 3-4


for (let i = 1; i <= job_array; i++) {
    // job index = 1-6

    

    console.log(`job_index: ${i}, row: ${row}, col: ${col}`);
}