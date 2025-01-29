// Removed duplicates from ACCRE output scraped text file.

const fs =  require('fs');
const outText = fs.readFileSync('./output/lig_adj_comp.txt', 'utf8');

const lines = outText.split('\n');
const uniqueLigands = [];

let realOut = ``;
lines.forEach(line => {
    const ligand1 = line.split(" ").slice(0, 2).join(" ");
    if (!uniqueLigands.includes(ligand1)) {
        uniqueLigands.push(ligand1);
        realOut += `${line}\n`;
    }
});

fs.writeFileSync('./output/lig_adj_comp_no_dupes.txt', realOut);
console.log(uniqueLigands.length);