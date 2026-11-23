struct Function {
    coefficients: Vec<f64>
}

impl Function {
    fn new() -> Self {
        Self {
            coefficients: Vec::new()
        }
    }

    fn parse(s: &str) -> Result<Function, String> {
        let mut monoms: Vec<(f64, usize)> = vec![];
        let mut curr_monom: (f64, usize) = (0f64, 0usize);
        let mut monoms_str: Vec<&str> = s.split_whitespace().collect();

        if monoms_str[0] != "+" && monoms_str[0] != "-" {
            monoms_str.insert(0usize, "+");
        }

        let mut max_pow: usize = 0;

        for i in 0..(monoms_str.len() / 2) {
            let sgn: f64 = if monoms_str[i * 2] == "-" { -1f64 } else { 1f64 };
            let monom_vec: Vec<&str> = monoms_str[i * 2 + 1].split("x").collect();

            match monom_vec[0].parse::<f64>() {
                Ok(num) => curr_monom.0 = num * sgn,
                Err(_) => {
                    if monom_vec[0] == "" { curr_monom.0 = sgn; }
                    else { return Err(String::from("Incorrect input!")); }
                }
            }
            
            if monom_vec.len() == 1 {
                curr_monom.1 = 0;
            }
            else {
                match monom_vec[1].parse::<usize>() {
                    Ok(num) => curr_monom.1 = num,
                    Err(_) => {
                        if monom_vec[1] == "" { curr_monom.1 = 1; }
                        else { return Err(String::from("Incorrect input!")); }
                    }
                }
            }

            max_pow = usize::max(max_pow, curr_monom.1);
            monoms.push(curr_monom);
        }

        let mut res: Function = Function::new();

        res.coefficients.resize(max_pow + 1, 0f64);

        for monom in monoms {
            res.coefficients[monom.1] += monom.0;
        }
        
        return Ok(res);
    }

    fn get_derivative(&self) -> Function {
        let mut derivative: Function = Function::new();
        derivative.coefficients.resize(self.coefficients.len() - 1, 0f64);

        for i in 0..(self.coefficients.len() - 1) {
            derivative.coefficients[i] = self.coefficients[i + 1] * ((i as f64) + 1f64);
        }

        derivative
    }

    fn degree(&self) -> i32 {
        self.coefficients.len() as i32 - 1
    }
    
    fn call(&self, x: f64) -> f64 {
        let mut y: f64 = 0f64;
        for i in 0..self.coefficients.len() {
            y += self.coefficients[i] * x.powi(i as i32);
        }
        y
    }

    fn get_root_binary(&self, mut a: f64, mut b: f64) -> Option<f64> {        
        if a > b { std::mem::swap(&mut a, &mut b); }
        if self.call(a) * self.call(b) > 0f64 || self.call(b) == 0f64 { return None; }
        if self.call(a) == 0f64 { return Some(a); }
        
        if self.call(a) > self.call(b) { std::mem::swap(&mut a, &mut b); }
        while (a - b).abs() > 1e-9f64 {
            let avg: f64 = (a + b) / 2f64;
            if self.call(avg) < 0f64 { a = avg; }
            else { b = avg }
        }
        
        let root: f64 = (a + b) / 2f64;
        if self.call(a) == 0f64 { return Some(a); }
        if self.call(b) == 0f64 { return Some(b); }
        for i in 0..4 {
            let pow: f64 = 10f64.powi(i);
            let rounded_root: f64 = (root * pow).round() / pow;
            if self.call((root * pow).round() / pow) == 0f64 && (rounded_root - a) * (rounded_root - b) < 0f64 { return Some(rounded_root); }
        }
        Some(root)
    }

    fn solve(&self) -> Vec<f64> {
        match self.degree() {
            1 => {
                let a: f64 = self.coefficients[1];
                let b: f64 = self.coefficients[0];
                if a == 0f64 {
                    return vec![];
                }
                return vec![-b / a];
            },
            2 => {
                let a: f64 = self.coefficients[2];
                let b: f64 = self.coefficients[1];
                let c: f64 = self.coefficients[0];

                let d: f64 = b.powi(2) - 4f64 * a * c;
                if d < 0f64 {
                    return vec![];
                }
                else if d == 0f64 {
                    return vec![-b / (2f64 * a)];
                }
                else {
                    if a > 0f64 {
                        vec![
                            (-b - d.sqrt()) / (2f64 * a),
                            (-b + d.sqrt()) / (2f64 * a),
                        ]
                    }
                    else {
                        vec![
                            (-b + d.sqrt()) / (2f64 * a),
                            (-b - d.sqrt()) / (2f64 * a),
                        ]
                    }
                }
            },
            _ => {
                let derivative: Function = self.get_derivative();
                let static_points: Vec<f64> = derivative.solve();

                let mut roots: Vec<f64> = Vec::new();

                if static_points.len() == 0 {
                    let mut border: f64 = 1f64;
                    while self.call(border) * self.call(-border) > 0f64 { border *= 2f64; }
                    match self.get_root_binary(-border, border) {
                        Some(root) => roots.push(root),
                        None => {}
                    }
                    return roots;
                }
                
                let sp_first: f64 = static_points[0];
                if self.call(sp_first) * (self.call(sp_first - 2f64 * sp_first.abs() - 1f64) - self.call(sp_first)) <= 0f64 {
                    let mut border: f64 = sp_first - 2f64 * sp_first.abs() - 1f64;
                    while self.call(sp_first) * self.call(border) > 0f64 { border *= 2f64; }
                    match self.get_root_binary(border, sp_first) {
                        Some(root) => roots.push(root),
                        None => {}
                    }
                }

                for i in 0..(static_points.len() - 1) {
                    match self.get_root_binary(static_points[i], static_points[i + 1]) {
                        Some(root) => roots.push(root),
                        None => {}
                    }
                }
                
                let sp_last: f64 = static_points[static_points.len() - 1];
                if self.call(sp_last) * (self.call(sp_last + 2f64 * sp_last.abs() + 1f64) - self.call(sp_last)) <= 0f64 {
                    let mut border: f64 = sp_last + 2f64 * sp_last.abs() + 1f64;
                    while self.call(sp_last) * self.call(border) > 0f64 { border *= 2f64; }
                    if self.call(border) == 0f64 { roots.push(border); }
                    else {
                        match self.get_root_binary(sp_last, border) {
                            Some(root) => roots.push(root),
                            None => {}
                        }
                    }
                }

                roots
            }
        }
    }
}

fn main() {
    let mut input: String = String::new();
    let mut function: Function;

    loop {
        input.clear();

        match std::io::stdin().read_line(&mut input) {
            Ok(_) => {},
            Err(e) => {
                println!("Terminal read error: {}\n", e);
                continue;
            }
        }

        match Function::parse(input.trim_matches(['\r', '\n'])) {
            Ok(f) => function = f,
            Err(e) => {
                println!("{}\n", e);
                continue;
            }
        }

        println!();

        let roots: Vec<f64> = function.solve();
        match roots.len() {
            0 => println!("This equation has no roots"),
            1 => println!("X = {}", roots[0]),
            _ => {
                for i in 0..roots.len() {
                    println!("X{} = {}", i + 1, roots[i]);
                }
            }
        }
        println!();
    }
}