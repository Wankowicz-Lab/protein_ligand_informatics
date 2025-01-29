// Scrape ACCRE output from ./output/out2 in case text file entries missing

const fs = require('fs');

const outputFiles = fs.readdirSync('./output/out2');
let newComp = ``;

function extractComparison(file) {
    const loaded = fs.readFileSync(`./output/out2/${file}`, 'utf8');
    const lines = loaded.split('\n');
    const comparison = lines.filter(line => line.includes('lig_adj'));
    let comp0 = comparison[0];
    if (comp0) {
        newComp += `${comp0}\n`;
    }
}

outputFiles.forEach(file => {
    extractComparison(file);
    console.log("running " + file)
});
console.log("done")
fs.writeFileSync('./output/lig_adj_comp.txt', newComp);