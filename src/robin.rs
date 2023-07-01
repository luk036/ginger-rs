#[derive(PartialEq, Eq, Clone, Debug, Default)]
pub struct SlNode {
    next: Option<Box<SlNode>>,
    data: usize,
}

pub struct RobinIterator<'a> {
    cur: &'a mut SlNode,
    stop: &'a SlNode,
}

impl<'a> Iterator for RobinIterator<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        self.cur = self.cur.next.as_mut().unwrap();
        if self.cur as *const SlNode != self.stop as *const SlNode {
            Some(self.cur.data)
        } else {
            None
        }
    }
}

pub struct Robin {
    cycle: Vec<SlNode>,
}

impl Robin {
    pub fn new(num_parts: usize) -> Self {
        let mut cycle = Vec::new();
        for k in 0..num_parts {
            cycle.push(SlNode {
                next: None,
                data: k,
            });
        }
        let mut sl2 = &mut cycle[num_parts - 1];
        for sl1 in &mut cycle {
            sl2.next = Some(Box::new(sl1.clone()));
            sl2 = sl2.next.as_mut().unwrap();
        }
        Self { cycle }
    }

    pub fn exclude(&mut self, from_part: usize) -> RobinIterator {
        RobinIterator {
            cur: &mut self.cycle[from_part],
            stop: &self.cycle[from_part],
        }
    }
}

// fn main() {
//     let mut r = Robin::new(5);
//     for k in r.exclude(3) {
//         println!("{}", k);
//     }
// }
