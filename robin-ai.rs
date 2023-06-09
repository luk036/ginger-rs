struct SlNode {
    next: Option<Box<SlNode>>,
    data: i32,
}

struct RobinIterator<'a> {
    cur: &'a mut SlNode,
    stop: &'a SlNode,
}

impl<'a> Iterator for RobinIterator<'a> {
    type Item = i32;

    fn next(&mut self) -> Option<Self::Item> {
        self.cur = self.cur.next.as_mut().unwrap();
        if self.cur as *const SlNode != self.stop as *const SlNode {
            Some(self.cur.data)
        } else {
            None
        }
    }
}

struct Robin {
    cycle: Vec<SlNode>,
}

impl Robin {
    fn new(num_parts: i32) -> Self {
        let mut cycle = Vec::new();
        for k in 0..num_parts {
            cycle.push(SlNode {
                next: None,
                data: k,
            });
        }
        let mut sl2 = &mut cycle[num_parts as usize - 1];
        for sl1 in &mut cycle {
            sl2.next = Some(Box::new(sl1.clone()));
            sl2 = sl2.next.as_mut().unwrap();
        }
        Self { cycle }
    }

    fn exclude(&mut self, from_part: usize) -> RobinIterator {
        RobinIterator {
            cur: &mut self.cycle[from_part],
            stop: &self.cycle[from_part],
        }
    }
}

fn main() {
    let mut r = Robin::new(5);
    for k in r.exclude(3) {
        println!("{}", k);
    }
}
