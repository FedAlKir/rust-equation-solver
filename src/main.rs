struct Function {
    coefficients: Vec<f64>
}

impl Function {
    fn new() -> Function {
        let new_f: Function = Function {
            coefficients: Vec::new()
        };
        new_f
    }

    fn parse(s: &String) -> Result<Function, String> {
        let mut monoms: Vec<(f64, usize)> = vec![];
        let mut curr_monom: (f64, usize) = (0f64, 0usize);
        let mut monoms_str: Vec<&str> = s.split(" ").collect();

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
                    if monom_vec[0] == "" { curr_monom.0 = 1f64; }
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
            res.coefficients[monom.1] = monom.0;
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

    fn get_pow(&self) -> i32 {
        return self.coefficients.len() as i32 - 1;
    }
    
    fn call(&self, x: f64) -> f64 {
        let mut y: f64 = 0f64;
        for i in 0..self.coefficients.len() {
            y += self.coefficients[i] * x.powi(i as i32);
        }
        y
    }

    fn get_root_binary(&self, mut a: f64, mut b: f64) -> f64 {
        if self.call(a) > self.call(b) { std::mem::swap(&mut a, &mut b); }
        while (a - b).abs() > 1e-9f64 {
            let avg: f64 = (a + b) / 2f64;
            if self.call(avg) < 0f64 { a = avg; }
            else { b = avg }
        }
        
        if self.call(a) == 0f64 { a }
        else if self.call(b) == 0f64 { b }
        else {
            let root: f64 = (a + b) / 2f64;
            for i in 0..4 {
                let pow: f64 = 10f64.powi(i);
                if self.call((root * pow).round() / pow) == 0f64 { return (root * pow).round() / pow; }
            }
            root
        }
    }

    fn solve(&self) -> Vec<f64> {
        match self.get_pow() {
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
                    return vec![
                        (-b - d.sqrt()) / (2f64 * a),
                        (-b + d.sqrt()) / (2f64 * a),
                    ]
                }
            },
            _ => {
                let derivative: Function = self.get_derivative();
                let static_points: Vec<f64> = derivative.solve();

                let mut roots: Vec<f64> = Vec::new();

                if static_points.len() == 0 {
                    let mut border: f64 = 1f64;
                    while self.call(border) * self.call(-border) > 0f64 { border *= 2f64; }
                    roots.push(self.get_root_binary(-border, border));
                }
                else {
                    let sp_first: f64 = static_points[0];
                    if self.call(sp_first) == 0f64 {
                        roots.push(sp_first);
                    }
                    else if self.call(sp_first) * (self.call(sp_first - 2f64 * sp_first.abs()) - self.call(sp_first)) < 0f64 {
                        let mut border: f64 = sp_first - 2f64 * sp_first.abs();
                        while self.call(sp_first) * self.call(border) > 0f64 { border *= 2f64}
                        roots.push(self.get_root_binary(border, sp_first));
                    }

                    for i in 0..(static_points.len() - 1) {
                        if self.call(static_points[i]) * self.call(static_points[i + 1]) > 0f64 {
                            continue;
                        }

                        roots.push(self.get_root_binary(static_points[i], static_points[i + 1]));
                    }
                    
                    if static_points.len() > 1 {
                        let sp_last: f64 = static_points[static_points.len() - 1];
                        if self.call(sp_last) == 0f64 {
                            roots.push(sp_last);
                        }
                        if self.call(sp_last) * (self.call(sp_last + 2f64 * sp_last.abs()) - self.call(sp_last)) < 0f64 {
                            let mut border: f64 = sp_last + 2f64 * sp_last.abs();
                            while self.call(sp_last) * self.call(border) > 0f64 { border *= 2f64}
                            roots.push(self.get_root_binary(sp_last, border));
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
        match std::io::stdin().read_line(&mut input) {
            Ok(_) => {},
            Err(e) => println!("Terminal read error: {}", e)
        }

        match Function::parse(&String::from(input.trim_matches(['\r', '\n']))) {
            Ok(f) => function = f,
            Err(e) => {
                println!("{}", e);
                continue;
            }
        }

        input.clear();

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
    }
}